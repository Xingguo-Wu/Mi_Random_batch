/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/vrrbl,PairLJCutVRRBL);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_VRRBL_H
#define LMP_PAIR_LJ_CUT_VRRBL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutVRRBL : public Pair {
 public:
  PairLJCutVRRBL(class LAMMPS *);
  ~PairLJCutVRRBL() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void born_matrix(int, int, int, int, double, double, double, double &, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double cut_global;
  double **cut;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;

  double core_size;
  int shell_samples;
  unsigned int random_seed_base;
  unsigned int random_state;

  int maxatom_cov;
  double **vrrbl_cov;

  void allocate();
  void allocate_covariance(int);
  unsigned int next_random();
};

}    // namespace LAMMPS_NS

#endif
#endif
