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

#include "pair_lj_cut_mirbl_lag.h"

#include <mpi.h>

#include <cmath>
#include <vector>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void PairLJCutMiRBLLag::settings(int narg, char **arg)
{
  if (narg != 7 && narg != 8)
    error->all(FLERR,
               "Illegal pair_style lj/cut/mirbl/lag command: expected cutoff core_size shell_samples rho alpha eps S0 [seed]");

  force_s0 = utils::numeric(FLERR, arg[6], false, lmp);
  if (force_s0 <= 0.0) error->all(FLERR, "Mi-RBL lag S0 must be > 0");

  std::vector<char *> base_args;
  base_args.reserve(7);
  for (int i = 0; i < 6; i++) base_args.push_back(arg[i]);
  if (narg == 8) base_args.push_back(arg[7]);
  PairLJCutMiRBL::settings(static_cast<int>(base_args.size()), base_args.data());
}
void PairLJCutMiRBLLag::compute(int eflag, int vflag)
{
  int i, ii, j, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  double rsq, r2inv, r6inv, forcelj, factor_lj;
  int *ilist, *jlist, *num_core_neigh, **first_core_neigh, *num_shell_neigh, **first_shell_neigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *mirbl_state = (alpha != 0.0) ? atom->dvector[mirbl_state_index] : nullptr;
  double *mirbl_gain = (mirbl_gain_index >= 0) ? atom->dvector[mirbl_gain_index] : nullptr;
  int *type = atom->type;
  double *special_lj = force->special_lj;

  inum = list->inum;
  ilist = list->ilist;
  num_core_neigh = list->num_core_neigh;
  first_core_neigh = list->first_core_neigh;
  num_shell_neigh = list->num_shell_neigh;
  first_shell_neigh = list->first_shell_neigh;

  std::vector<int> shell_choice(shell_samples);
  const double ema_weight = (alpha != 0.0) ? 1.0 - rho : 0.0;
  const double whitening_power = 0.5 * alpha;
  const bool has_mirbl_history = (alpha != 0.0 && mirbl_updates > 0);
  double safe_bias_correction = 1.0;
  if (has_mirbl_history) {
    const double bias_correction = 1.0 - pow(rho, static_cast<double>(mirbl_updates));
    safe_bias_correction = (bias_correction > 1.0e-16) ? bias_correction : 1.0e-16;
  }

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

    // Shell neighbors are sampled, unbiasedly rescaled, and isotropically normalized.
    jlist = first_shell_neigh[i];
    jnum = num_shell_neigh[i];
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

    double shell_force[3] = {0.0, 0.0, 0.0};
    for (int isample = 0; isample < nsample; isample++) {
      j = shell_choice[isample];
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

        shell_force[0] += delx * fpair;
        shell_force[1] += dely * fpair;
        shell_force[2] += delz * fpair;
      }
    }

    const double unbiased_scale = (nsample > 0) ? static_cast<double>(jnum) / nsample : 0.0;
    shell_force[0] *= unbiased_scale;
    shell_force[1] *= unbiased_scale;
    shell_force[2] *= unbiased_scale;

    double mirbl_force[3] = {shell_force[0], shell_force[1], shell_force[2]};
    double gain = 1.0;
    if (alpha != 0.0) {
      const double force2_over3 =
          (shell_force[0] * shell_force[0] + shell_force[1] * shell_force[1] +
           shell_force[2] * shell_force[2]) /
          3.0;

      if (has_mirbl_history) {
        const double corrected_state = mirbl_state[i] / safe_bias_correction;
        gain = pow((corrected_state + eps) / force_s0, -whitening_power);
      }
      mirbl_state[i] = rho * mirbl_state[i] + ema_weight * force2_over3;

      mirbl_force[0] = gain * shell_force[0];
      mirbl_force[1] = gain * shell_force[1];
      mirbl_force[2] = gain * shell_force[2];
    }
    if (mirbl_gain) mirbl_gain[i] = gain;

    f[i][0] += mirbl_force[0];
    f[i][1] += mirbl_force[1];
    f[i][2] += mirbl_force[2];

    sampled_force_sum[0] += mirbl_force[0];
    sampled_force_sum[1] += mirbl_force[1];
    sampled_force_sum[2] += mirbl_force[2];
  }

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

  if (alpha != 0.0) mirbl_updates++;
}

/* ---------------------------------------------------------------------- */

void PairLJCutMiRBLLag::init_style()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Mi-RBL lag pair style does not yet support run_style respa");

  if (alpha != 0.0) {
    int flag = -1;
    int cols = -1;
    mirbl_state_index = atom->find_custom("mirbl_state", flag, cols);
    if (mirbl_state_index < 0 || flag != 1 || cols != 0)
      error->all(FLERR,
                 "Mi-RBL lag pair style requires 'fix <id> all property/atom d_mirbl_state'");
  }

  int flag, cols;
  mirbl_gain_index = atom->find_custom("mirbl_gain", flag, cols);

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_CORE_SHELL_ONLY);
}
