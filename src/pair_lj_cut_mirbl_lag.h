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
   Contributing authors: Xingguo Wu (SJTU)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/mirbl/lag,PairLJCutMiRBLLag);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_MIRBL_LAG_H
#define LMP_PAIR_LJ_CUT_MIRBL_LAG_H

#include "pair_lj_cut_mirbl.h"

namespace LAMMPS_NS {

class PairLJCutMiRBLLag : public PairLJCutMiRBL {
 public:
  PairLJCutMiRBLLag(class LAMMPS *lmp) : PairLJCutMiRBL(lmp), force_s0(1.0), mirbl_gain_index(-1) {}
  void compute(int, int) override;
  void settings(int, char **) override;
  void init_style() override;

 protected:
  double force_s0;
  int mirbl_gain_index;
};

}    // namespace LAMMPS_NS

#endif
#endif
