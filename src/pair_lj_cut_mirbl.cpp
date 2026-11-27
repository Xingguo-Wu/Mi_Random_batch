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
   Contributing authors: Xingguo Wu (SJTU)
------------------------------------------------------------------------- */

#include "pair_lj_cut_mirbl.h"

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

/* ---------------------------------------------------------------------- */

PairLJCutMiRBL::PairLJCutMiRBL(LAMMPS *lmp) : Pair(lmp)
{
  born_matrix_enable = 1;
  writedata = 1;

  cut_global = 0.0;
  core_size = 0.0;
  shell_samples = 0;
  rho = 0.0;
  alpha = 0.0;
  eps = 1.0e-8;
  random_seed_base = 12345U;
  random_state = random_seed_base + static_cast<unsigned int>(comm->me);
  mirbl_state_index = -1;
  mirbl_updates = 0;
}

/* ---------------------------------------------------------------------- */

PairLJCutMiRBL::~PairLJCutMiRBL()
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
}

/* ---------------------------------------------------------------------- */

unsigned int PairLJCutMiRBL::next_random()
{
  random_state = random_state * 1103515245U + 12345U;
  return random_state;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutMiRBL::allocate()
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

void PairLJCutMiRBL::settings(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  core_size = utils::numeric(FLERR, arg[1], false, lmp);
  shell_samples = utils::inumeric(FLERR, arg[2], false, lmp);
  rho = utils::numeric(FLERR, arg[3], false, lmp);
  alpha = utils::numeric(FLERR, arg[4], false, lmp);
  eps = utils::numeric(FLERR, arg[5], false, lmp);
  if (narg == 7) {
    const int parsed_seed = utils::inumeric(FLERR, arg[6], false, lmp);
    if (parsed_seed < 0) error->all(FLERR, "Mi-RBL random seed must be >= 0");
    random_seed_base = static_cast<unsigned int>(parsed_seed);
  }
  random_state = random_seed_base + static_cast<unsigned int>(comm->me);
  mirbl_updates = 0;

  if (core_size < 0.0) error->all(FLERR, "Mi-RBL core_size must be >= 0");
  if (shell_samples <= 0) error->all(FLERR, "Mi-RBL shell sample count must be > 0");
  if (rho < 0.0 || rho >= 1.0) error->all(FLERR, "Mi-RBL rho must satisfy 0 <= rho < 1");
  if (alpha < 0.0 || alpha > 1.0) error->all(FLERR, "Mi-RBL alpha must satisfy 0 <= alpha <= 1");
  if (eps <= 0.0) error->all(FLERR, "Mi-RBL eps must be > 0");

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  force->pair->cut_core = core_size;

  if (force->newton_pair == 1)
    error->all(FLERR, "Mi-RBL pair styles require newton pair off");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutMiRBL::coeff(int narg, char **arg)
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
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutMiRBL::init_one(int i, int j)
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

void PairLJCutMiRBL::write_restart(FILE *fp)
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

void PairLJCutMiRBL::read_restart(FILE *fp)
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

void PairLJCutMiRBL::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
  fwrite(&core_size, sizeof(double), 1, fp);
  fwrite(&shell_samples, sizeof(int), 1, fp);
  fwrite(&rho, sizeof(double), 1, fp);
  fwrite(&alpha, sizeof(double), 1, fp);
  fwrite(&eps, sizeof(double), 1, fp);
  fwrite(&random_seed_base, sizeof(unsigned int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutMiRBL::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &core_size, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &shell_samples, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &rho, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &alpha, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &eps, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &random_seed_base, sizeof(unsigned int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&core_size, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&shell_samples, 1, MPI_INT, 0, world);
  MPI_Bcast(&rho, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&random_seed_base, sizeof(unsigned int), MPI_BYTE, 0, world);

  random_state = random_seed_base + static_cast<unsigned int>(comm->me);
  mirbl_updates = 0;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCutMiRBL::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g\n", i, epsilon[i][i], sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCutMiRBL::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJCutMiRBL::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
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

void PairLJCutMiRBL::born_matrix(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
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

void *PairLJCutMiRBL::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  return nullptr;
}
