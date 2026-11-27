/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: VR-RBL method: Qi Zhou (SJTU)
                         VR-RBL implementation: Xingguo Wu (SJTU)
------------------------------------------------------------------------- */

#include "pair_lj_cut_vrrbl.h"

#include <mpi.h>

#include <cmath>
#include <cstring>
#include <vector>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

namespace {

enum { COV_XX = 0, COV_XY = 1, COV_XZ = 2, COV_YY = 3, COV_YZ = 4, COV_ZZ = 5 };

}

/* ---------------------------------------------------------------------- */

PairLJCutVRRBL::PairLJCutVRRBL(LAMMPS *lmp) : Pair(lmp)
{
  born_matrix_enable = 1;
  writedata = 1;

  cut_global = 0.0;
  core_size = 0.0;
  shell_samples = 0;
  random_seed_base = 12345U;
  random_state = random_seed_base + static_cast<unsigned int>(comm->me);

  maxatom_cov = 0;
  vrrbl_cov = nullptr;
}

/* ---------------------------------------------------------------------- */

PairLJCutVRRBL::~PairLJCutVRRBL()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }

  memory->destroy(vrrbl_cov);
}

/* ---------------------------------------------------------------------- */

void PairLJCutVRRBL::allocate_covariance(int nmax)
{
  if (nmax <= maxatom_cov) return;
  memory->destroy(vrrbl_cov);
  maxatom_cov = nmax;
  memory->create(vrrbl_cov, maxatom_cov, 6, "pair:vrrbl_cov");
}

/* ---------------------------------------------------------------------- */

unsigned int PairLJCutVRRBL::next_random()
{
  random_state = random_state * 1103515245U + 12345U;
  return random_state;
}

/* ----------------------------------------------------------------------
   compute pure-LJ pairwise forces using exact core neighbors plus a
   sampled shell estimator. Alongside the sampled shell force, estimate
   the instantaneous covariance needed by fix vr/langevin.
------------------------------------------------------------------------- */

void PairLJCutVRRBL::compute(int eflag, int vflag)
{
  int i, ii, j, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  double rsq, r2inv, r6inv, forcelj, factor_lj;
  int *ilist, *jlist, *num_core_neigh, **first_core_neigh, *num_shell_neigh, **first_shell_neigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  double *special_lj = force->special_lj;

  inum = list->inum;
  ilist = list->ilist;
  num_core_neigh = list->num_core_neigh;
  first_core_neigh = list->first_core_neigh;
  num_shell_neigh = list->num_shell_neigh;
  first_shell_neigh = list->first_shell_neigh;

  allocate_covariance(atom->nmax);
  for (i = 0; i < atom->nlocal; i++)
    for (int cov = 0; cov < 6; cov++) vrrbl_cov[i][cov] = 0.0;

  std::vector<int> shell_choice(shell_samples);
  double sampled_force_sum[3] = {0.0, 0.0, 0.0};

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    // Core neighbors are accumulated exactly.
    jlist = first_core_neigh[i];
    jnum = num_core_neigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;
        r6inv = r2inv * r2inv * r2inv;
        forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
        fpair = factor_lj * forcelj * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
      }
    }

    // Shell neighbors are sampled without replacement, matching the
    // RBL/Mi-RBL shell sampling convention.
    jlist = first_shell_neigh[i];
    jnum = num_shell_neigh[i];
    if (jnum == 0) continue;
    const int nsample = (jnum < shell_samples) ? jnum : shell_samples;

    int remaining = nsample;
    int picked = 0;
    for (int ishell = 0; ishell < jnum; ishell++) {
      if (remaining == 0) break;
      if (next_random() % static_cast<unsigned int>(jnum - ishell) <
          static_cast<unsigned int>(remaining)) {
        shell_choice[picked++] = jlist[ishell];
        remaining--;
      }
    }

    double shell_force_sum[3] = {0.0, 0.0, 0.0};
    double shell_second[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Previous with-replacement sampler kept for comparison/recovery:
    // for (int isample = 0; isample < shell_samples; isample++) {
    //   j = jlist[next_random() % static_cast<unsigned int>(jnum)];
    //   factor_lj = special_lj[sbmask(j)];
    //   j &= NEIGHMASK;
    //
    //   delx = xtmp - x[j][0];
    //   dely = ytmp - x[j][1];
    //   delz = ztmp - x[j][2];
    //   rsq = delx * delx + dely * dely + delz * delz;
    //   jtype = type[j];
    //
    //   if (rsq >= cutsq[itype][jtype]) continue;
    //
    //   r2inv = 1.0 / rsq;
    //   r6inv = r2inv * r2inv * r2inv;
    //   forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
    //   fpair = factor_lj * forcelj * r2inv;
    //
    //   const double fx = delx * fpair;
    //   const double fy = dely * fpair;
    //   const double fz = delz * fpair;
    //
    //   shell_force_sum[0] += fx;
    //   shell_force_sum[1] += fy;
    //   shell_force_sum[2] += fz;
    //
    //   shell_second[COV_XX] += fx * fx;
    //   shell_second[COV_XY] += fx * fy;
    //   shell_second[COV_XZ] += fx * fz;
    //   shell_second[COV_YY] += fy * fy;
    //   shell_second[COV_YZ] += fy * fz;
    //   shell_second[COV_ZZ] += fz * fz;
    // }

    for (int isample = 0; isample < nsample; isample++) {
      j = shell_choice[isample];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      r2inv = 1.0 / rsq;
      r6inv = r2inv * r2inv * r2inv;
      forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
      fpair = factor_lj * forcelj * r2inv;

      const double fx = delx * fpair;
      const double fy = dely * fpair;
      const double fz = delz * fpair;

      shell_force_sum[0] += fx;
      shell_force_sum[1] += fy;
      shell_force_sum[2] += fz;

      shell_second[COV_XX] += fx * fx;
      shell_second[COV_XY] += fx * fy;
      shell_second[COV_XZ] += fx * fz;
      shell_second[COV_YY] += fy * fy;
      shell_second[COV_YZ] += fy * fz;
      shell_second[COV_ZZ] += fz * fz;
    }

    const double nsample_double = static_cast<double>(nsample);
    const double scale = static_cast<double>(jnum) / nsample_double;
    // Paper-faithful VR-RBL covariance prefactor from Eq. (2.8), kept for
    // main-text baseline runs:
    // const double cov_scale =
    //     (nsample < jnum) ?
    //     static_cast<double>(jnum) * static_cast<double>(jnum) /
    //         (nsample_double * nsample_double) :
    //     0.0;

    // Finite-population correction for the current no-replacement ablation.
    const double cov_scale =
        (nsample > 1 && nsample < jnum) ?
        static_cast<double>(jnum) * static_cast<double>(jnum - nsample) /
            (nsample_double * static_cast<double>(nsample - 1)) :
        0.0;
    const double scaled_shell_force[3] = {
        shell_force_sum[0] * scale, shell_force_sum[1] * scale, shell_force_sum[2] * scale};

    f[i][0] += scaled_shell_force[0];
    f[i][1] += scaled_shell_force[1];
    f[i][2] += scaled_shell_force[2];

    sampled_force_sum[0] += scaled_shell_force[0];
    sampled_force_sum[1] += scaled_shell_force[1];
    sampled_force_sum[2] += scaled_shell_force[2];

    vrrbl_cov[i][COV_XX] =
        cov_scale * (shell_second[COV_XX] - shell_force_sum[0] * shell_force_sum[0] / nsample);
    vrrbl_cov[i][COV_XY] =
        cov_scale * (shell_second[COV_XY] - shell_force_sum[0] * shell_force_sum[1] / nsample);
    vrrbl_cov[i][COV_XZ] =
        cov_scale * (shell_second[COV_XZ] - shell_force_sum[0] * shell_force_sum[2] / nsample);
    vrrbl_cov[i][COV_YY] =
        cov_scale * (shell_second[COV_YY] - shell_force_sum[1] * shell_force_sum[1] / nsample);
    vrrbl_cov[i][COV_YZ] =
        cov_scale * (shell_second[COV_YZ] - shell_force_sum[1] * shell_force_sum[2] / nsample);
    vrrbl_cov[i][COV_ZZ] =
        cov_scale * (shell_second[COV_ZZ] - shell_force_sum[2] * shell_force_sum[2] / nsample);
  }

  // Match the existing RBL/Mi-RBL convention: remove the mean sampled force.
  double sampled_force_sum_all[3] = {0.0, 0.0, 0.0};
  MPI_Allreduce(sampled_force_sum, sampled_force_sum_all, 3, MPI_DOUBLE, MPI_SUM, world);

  const double inv_natoms = 1.0 / static_cast<double>(atom->natoms);
  const double mean_sampled_force[3] = {sampled_force_sum_all[0] * inv_natoms,
                                        sampled_force_sum_all[1] * inv_natoms,
                                        sampled_force_sum_all[2] * inv_natoms};

  for (i = 0; i < atom->nlocal; i++) {
    f[i][0] -= mean_sampled_force[0];
    f[i][1] -= mean_sampled_force[1];
    f[i][2] -= mean_sampled_force[2];
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutVRRBL::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(cut, n, n, "pair:cut");
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(sigma, n, n, "pair:sigma");
  memory->create(lj1, n, n, "pair:lj1");
  memory->create(lj2, n, n, "pair:lj2");
  memory->create(lj3, n, n, "pair:lj3");
  memory->create(lj4, n, n, "pair:lj4");
  memory->create(offset, n, n, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutVRRBL::settings(int narg, char **arg)
{
  if (narg != 3 && narg != 4) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  core_size = utils::numeric(FLERR, arg[1], false, lmp);
  shell_samples = utils::inumeric(FLERR, arg[2], false, lmp);
  if (narg == 4) {
    const int parsed_seed = utils::inumeric(FLERR, arg[3], false, lmp);
    if (parsed_seed < 0) error->all(FLERR, "VRRBL random seed must be >= 0");
    random_seed_base = static_cast<unsigned int>(parsed_seed);
  }
  random_state = random_seed_base + static_cast<unsigned int>(comm->me);

  if (core_size < 0.0) error->all(FLERR, "VRRBL core_size must be >= 0");
  if (shell_samples <= 0) error->all(FLERR, "VRRBL shell sample count must be > 0");

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  force->pair->cut_core = core_size;

  if (force->newton_pair == 1)
    error->all(FLERR, "Pair style lj/cut/vrrbl requires newton pair off");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutVRRBL::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);

  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutVRRBL::init_style()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Pair style lj/cut/vrrbl does not yet support run_style respa");

  allocate_covariance(atom->nmax);
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_CORE_SHELL_ONLY);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutVRRBL::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
  } else
    offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2], all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    double sig2 = sigma[i][j] * sigma[i][j];
    double sig6 = sig2 * sig2 * sig2;
    double rc3 = cut[i][j] * cut[i][j] * cut[i][j];
    double rc6 = rc3 * rc3;
    double rc9 = rc3 * rc6;
    double prefactor = 8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 / (9.0 * rc9);
    etail_ij = prefactor * (sig6 - 3.0 * rc6);
    ptail_ij = 2.0 * prefactor * (2.0 * sig6 - 3.0 * rc6);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutVRRBL::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutVRRBL::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int me = comm->me;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutVRRBL::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
  fwrite(&core_size, sizeof(double), 1, fp);
  fwrite(&shell_samples, sizeof(int), 1, fp);
  fwrite(&random_seed_base, sizeof(unsigned int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutVRRBL::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &core_size, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &shell_samples, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &random_seed_base, sizeof(unsigned int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&core_size, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&shell_samples, 1, MPI_INT, 0, world);
  MPI_Bcast(&random_seed_base, sizeof(unsigned int), MPI_BYTE, 0, world);

  random_state = random_seed_base + static_cast<unsigned int>(comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCutVRRBL::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g\n", i, epsilon[i][i], sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCutVRRBL::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJCutVRRBL::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                              double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r2inv, r6inv, forcelj, philj;

  r2inv = 1.0 / rsq;
  r6inv = r2inv * r2inv * r2inv;
  forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
  fforce = factor_lj * forcelj * r2inv;

  philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
  return factor_lj * philj;
}

/* ---------------------------------------------------------------------- */

void PairLJCutVRRBL::born_matrix(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                                 double /*factor_coul*/, double factor_lj, double &dupair,
                                 double &du2pair)
{
  double rinv, r2inv, r6inv, du, du2;

  r2inv = 1.0 / rsq;
  rinv = sqrt(r2inv);
  r6inv = r2inv * r2inv * r2inv;

  du = r6inv * rinv * (lj2[itype][jtype] - lj1[itype][jtype] * r6inv);
  du2 = r6inv * r2inv * (13 * lj1[itype][jtype] * r6inv - 7 * lj2[itype][jtype]);

  dupair = factor_lj * du;
  du2pair = factor_lj * du2;
}

/* ---------------------------------------------------------------------- */

void *PairLJCutVRRBL::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  if (strcmp(str, "vrrbl_cov") == 0) return (void *) vrrbl_cov;
  return nullptr;
}
