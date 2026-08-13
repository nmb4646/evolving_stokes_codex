#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_dir="${script_dir}/vendor/mmg"
build_dir="${script_dir}/build"
install_dir="${script_dir}/install"
mmg_tag="v5.8.0"
mmg_commit="4d8232c8aebfed877935d75d4d4a67e850962422"
matlab_bin="${MATLAB_BIN:-matlab}"

if [[ ! -d "${source_dir}/.git" ]]; then
    mkdir -p "$(dirname "${source_dir}")"
    git clone --depth 1 --branch "${mmg_tag}" \
        https://github.com/MmgTools/mmg.git "${source_dir}"
fi

actual_commit="$(git -C "${source_dir}" rev-parse HEAD)"
if [[ "${actual_commit}" != "${mmg_commit}" ]]; then
    echo "Expected MMG ${mmg_tag} at ${mmg_commit}, found ${actual_commit}." >&2
    echo "Remove ${source_dir} and rerun this script to restore the pinned source." >&2
    exit 1
fi

cmake -S "${source_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${install_dir}" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD=MMGS \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_STATIC_LIBS=ON \
    -DUSE_SCOTCH=OFF \
    -DBUILD_TESTING=OFF \
    -DTEST_LIBMMGS=OFF

cmake --build "${build_dir}" --parallel "$(nproc)"
cmake --install "${build_dir}"

"${matlab_bin}" -nodesktop -nosplash -batch \
    "mex('-R2018a','-I${install_dir}/include','${script_dir}/mmgs_remesh_mex.c','${install_dir}/lib/libmmgs.a','-lstdc++','-outdir','${script_dir}')"

echo "Built ${script_dir}/mmgs_remesh_mex.mexa64 from MMG ${mmg_tag}."
