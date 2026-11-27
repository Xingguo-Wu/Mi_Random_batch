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

#include "npair_full_bin_atomonly_coreshell.h"

#include "atom.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairFullBinAtomonlyCoreShell::NPairFullBinAtomonlyCoreShell(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   build only the RBL core/shell full lists for atom-only systems
------------------------------------------------------------------------- */

void NPairFullBinAtomonlyCoreShell::build(NeighList *list)
{
  int i, j, k, n_core, n_shell, itype, jtype, ibin;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *neigh_core_ptr;
  int *neigh_shell_ptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *num_core_neigh = list->num_core_neigh;
  int **first_core_neigh = list->first_core_neigh;
  int *num_shell_neigh = list->num_shell_neigh;
  int **first_shell_neigh = list->first_shell_neigh;
  MyPage<int> *ipage_core = list->ipage_core;
  MyPage<int> *ipage_shell = list->ipage_shell;

  int inum = 0;
  ipage_core->reset();
  ipage_shell->reset();

  for (i = 0; i < nlocal; i++) {
    n_core = 0;
    n_shell = 0;
    neigh_core_ptr = ipage_core->vget();
    neigh_shell_ptr = ipage_shell->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    const double core_size = cut_core_;
    const double shell_size = cut_shell_;

    ibin = atom2bin[i];

    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
        if (i == j) continue;

        jtype = type[j];
        if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        if (rsq <= cutneighsq[itype][jtype]) {
          if (rsq <= core_size * core_size)
            neigh_core_ptr[n_core++] = j;
          else if (rsq <= shell_size * shell_size)
            neigh_shell_ptr[n_shell++] = j;
        }
      }
    }

    ilist[inum++] = i;
    first_core_neigh[i] = neigh_core_ptr;
    num_core_neigh[i] = n_core;
    first_shell_neigh[i] = neigh_shell_ptr;
    num_shell_neigh[i] = n_shell;

    ipage_core->vgot(n_core);
    if (ipage_core->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

    ipage_shell->vgot(n_shell);
    if (ipage_shell->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
  list->gnum = 0;
}
