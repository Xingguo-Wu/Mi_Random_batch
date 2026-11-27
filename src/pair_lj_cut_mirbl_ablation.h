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
PairStyle(lj/cut/mirbl/ablation/global,PairLJCutMiRBLAblationGlobal);
PairStyle(lj/cut/mirbl/ablation/species,PairLJCutMiRBLAblationSpecies);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_MIRBL_ABLATION_H
#define LMP_PAIR_LJ_CUT_MIRBL_ABLATION_H

#include "pair_lj_cut_mirbl.h"

#include <vector>

namespace LAMMPS_NS {

class PairLJCutMiRBLAblation : public PairLJCutMiRBL {
 public:
  enum Mode { GLOBAL, SPECIES };

  PairLJCutMiRBLAblation(class LAMMPS *lmp, Mode mode) : PairLJCutMiRBL(lmp), mode(mode), force_s0(1.0), mirbl_gain_index(-1) {}
  void compute(int, int) override;
  void settings(int, char **) override;
  void init_style() override;

 protected:
  Mode mode;
  double force_s0;
  int mirbl_gain_index;
  std::vector<double> group_state;
};

class PairLJCutMiRBLAblationGlobal : public PairLJCutMiRBLAblation {
 public:
  PairLJCutMiRBLAblationGlobal(class LAMMPS *lmp) : PairLJCutMiRBLAblation(lmp, GLOBAL) {}
};

class PairLJCutMiRBLAblationSpecies : public PairLJCutMiRBLAblation {
 public:
  PairLJCutMiRBLAblationSpecies(class LAMMPS *lmp) : PairLJCutMiRBLAblation(lmp, SPECIES) {}
};

}    // namespace LAMMPS_NS

#endif
#endif
