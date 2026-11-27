#!/usr/bin/env bash
set -euo pipefail

: "${LAMMPS_BIN:?Set LAMMPS_BIN to the absolute path of the LAMMPS executable}"
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
run_dir=${1:-example_run}
mkdir -p "$run_dir"
cd "$run_dir"

"$LAMMPS_BIN" -in "$script_dir/in.binary_lj" -log log.lammps \
  -var ncells 5 -var equil_steps 10 -var prod_steps 20 \
  -var thermo_every 5 -var sample_every 5 -var average_count 1 \
  -var traj_every 5 -var gain_every 5
