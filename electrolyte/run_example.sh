#!/usr/bin/env bash
set -euo pipefail

: "${LAMMPS_BIN:?Set LAMMPS_BIN to the absolute path of the LAMMPS executable}"
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
run_dir=${1:-example_run}
mkdir -p "$run_dir"
cd "$run_dir"

"$LAMMPS_BIN" -in "$script_dir/in.electrolyte" -log log.lammps \
  -var n_cation 100 -var n_anion 300 -var lbox 12.8 \
  -var cut_coul 4.0 \
  -var equil_steps 10 -var prod_steps 20 -var thermo_every 5 \
  -var traj_every 5 -var gain_every 5
