#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
build_dir="$project_root/build"

mkdir -p "$build_dir"
cd "$build_dir"

cmake "$project_root/cmake" \
  -D BUILD_MPI=ON \
  -D PKG_KSPACE=ON

cmake --build "$build_dir"

echo "LAMMPS executable: $build_dir/lmp"
