# Full Table 1 reproduction (26 circuits, crash-safe per-circuit runs).
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)

$image = "qsyn-test-gcc"
if (-not (docker image inspect $image 2>$null)) {
    docker build -t $image -f docker/gcc-test.Dockerfile .
}

Write-Host "Building + running full Table 1 (may take 1-2 hours for large circuits)..." -ForegroundColor Cyan
docker run --rm `
    --entrypoint /bin/bash `
    --security-opt seccomp=unconfined `
    -e TABLE1_TIMEOUT=1200 `
    -v "${PWD}:/app/qsyn" `
    -v qsyn-build-cache:/tmp/qsyn-build `
    $image `
    -lc "cp /app/qsyn/docker/test-entrypoint.sh /tmp/e.sh && sed -i 's/\r$//' /tmp/e.sh /app/qsyn/scripts/reproduce_table1.sh && QSYN_RUN_INTEGRATION=0 QSYN_RUN_DEMO=0 bash /tmp/e.sh '[phasepoly]' >/dev/null 2>&1 && cmake --build /tmp/qsyn-build --target qsyn -j >/dev/null && bash /app/qsyn/scripts/reproduce_table1.sh"

if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Results: doc/table1_results.md and doc/table1_results.csv" -ForegroundColor Green
