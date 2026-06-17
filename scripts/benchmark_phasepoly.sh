#!/usr/bin/env bash
# Run PhasePoly integration demo inside Docker (from repo root).
set -euo pipefail
cd "$(dirname "$0")/.."
exec docker run --rm --entrypoint /bin/bash --security-opt seccomp=unconfined \
  -v "$(pwd):/app/qsyn" qsyn-test-gcc -lc \
  "find /app/qsyn -name '*.sh' -exec sed -i 's/\r$//' {} + 2>/dev/null; \
   sed -i 's/\r$//' /app/entrypoint.sh 2>/dev/null; \
   cmake -B /app/build -S /app/qsyn -DCMAKE_BUILD_TYPE=Release && \
   cmake --build /app/build --target qsyn --parallel \$(nproc) && \
   cd /app/qsyn && /app/build/qsyn -v examples/phasepoly_spidernest.dof"
