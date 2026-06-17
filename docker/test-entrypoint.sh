#! /usr/bin/env bash
# Build in /tmp (not the bind-mounted repo) to avoid duplicating artifacts on the host.
set -euo pipefail

BUILD_DIR="${QSYN_BUILD_DIR:-/tmp/qsyn-build}"
SRC="/app/qsyn"

# GCC 12 + libpopcnt AVX512 triggers -Werror=uninitialized in system headers.
export CXXFLAGS="${CXXFLAGS:-} -Wno-error=uninitialized"

cmake -B "${BUILD_DIR}" -S "${SRC}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --parallel "$(nproc)" --target unit-test qsyn

cd "${SRC}" || exit 1

if [[ $# -gt 0 ]]; then
  "${BUILD_DIR}/qsyn-unit-test" "$@"
else
  "${BUILD_DIR}/qsyn-unit-test"
fi

if [[ "${QSYN_RUN_INTEGRATION:-0}" == "1" ]]; then
  "${SRC}/scripts/RUN_TESTS" --qsyn "${BUILD_DIR}/qsyn"
fi

# Optional demo when QSYN_RUN_DEMO=1
if [[ "${QSYN_RUN_DEMO:-0}" == "1" ]]; then
  "${BUILD_DIR}/qsyn" -v "${SRC}/examples/phasepoly_spidernest.dof"
fi
