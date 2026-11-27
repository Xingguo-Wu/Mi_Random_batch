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

#ifdef KSPACE_CLASS

KSpaceStyle(mi/rbe/lag,MiRbe)

#else

#ifndef LMP_MI_RBE_H
#define LMP_MI_RBE_H

#include "kspace.h"
#include<iostream>
#include <fstream>
#include <sstream>

namespace LAMMPS_NS {
	
class MiRbe : public KSpace {
public:
  MiRbe(class LAMMPS *);
  virtual ~MiRbe();
  void settings(int, char**) override;
  void init() override;
  void setup() override;
  void compute(int, int) override;
  double memory_usage() override;
  void slabcorr();

protected:
  int triclinic;
  double alpha, alpha_sqrt;
  int P;                    
  double S;
  double pi,pi_sqrt;
  double volume;
  std::streamoff pointer; 
  std::ifstream infile;
  int lx,ly,lz;
  int Step;
  int me;
  int (*K_All)[3];
  int RankID;
  double S1_X, S1_Y, S1_Z;
  double S_X, S_Y, S_Z;
  double xprd0, yprd0, zprd0;
  double S0;
  int** K_Sample;
  double* TimeSet;
  double mirbe_rho, mirbe_alpha, mirbe_eps;
  double force_s0;
  unsigned int random_seed_base;
  bigint mirbe_updates;
  int mirbe_state_index;
  int mirbe_gain_index;
  int write_timing;
};
}
#endif
#endif

