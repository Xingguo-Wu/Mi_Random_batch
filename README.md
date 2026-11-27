# Moment-informed force rescaling for random-batch Langevin dynamics

This research repository is intended for reproducing the paper results and studying the methods. It contains a LAMMPS 15 Jun 2023 based implementation of Mi-RBL and Mi-RBE, VR-RBL comparison code, simulation inputs, processed paper data, and plotting scripts.

## Repository layout

- `src/`, `cmake/`, `doc/`: LAMMPS source code, build files, and upstream documentation.
- `coexist_lj/`, `binary_lj/`, `electrolyte/`: inputs and short instructions for the three simulation systems.
- `scripts/compile.sh`: a general CMake build script with MPI and KSPACE enabled.
- `paper/data/`, `paper/plotting/`, `paper/postprocessing/`, `paper/figures/`: processed data, plotting scripts, postprocessing scripts, and generated figures.

## Getting started

With CMake, a C++ compiler, and MPI available, run `bash scripts/compile.sh` from the repository root. The executable is written to `build/lmp`. See the README in each simulation directory for its input and run instructions. Production numerical experiments are run on a cluster.

## Citation and license

Please cite the associated manuscript, *Moment-informed force rescaling for random-batch Langevin dynamics* (Xingguo Wu, Yangshuai Wang, and Zhenli Xu), and LAMMPS: [Thompson et al., Computer Physics Communications 271, 108171 (2022)](https://doi.org/10.1016/j.cpc.2021.108171). Citation metadata is in `CITATION.cff`; the original LAMMPS citation file is preserved as `CITATION.LAMMPS.cff`.

The LAMMPS-derived source is distributed under GPL v2; see `LICENSE`. Original LAMMPS information is preserved in `README.LAMMPS`, and bundled components retain their own license notices.
