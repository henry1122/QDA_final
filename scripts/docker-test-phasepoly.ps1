# Lightweight PhasePoly test via Docker (single build in /tmp, cached volume).
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)

$image = "qsyn-test-gcc"
if (-not (docker image inspect $image 2>$null)) {
    Write-Host "Building $image (one-time)..." -ForegroundColor Cyan
    docker build -t $image -f docker/gcc-test.Dockerfile .
}

Write-Host "Building + running PhasePoly tests and spidernest demo..." -ForegroundColor Cyan
docker run --rm `
    --entrypoint /bin/bash `
    --security-opt seccomp=unconfined `
    -e QSYN_RUN_INTEGRATION=0 `
    -e QSYN_RUN_DEMO=1 `
    -v "${PWD}:/app/qsyn:ro" `
    -v qsyn-build-cache:/tmp/qsyn-build `
    $image `
    -lc "sed -i 's/\r$//' /app/qsyn/docker/test-entrypoint.sh 2>/dev/null; bash /app/qsyn/docker/test-entrypoint.sh '[phasepoly]'"

if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Done." -ForegroundColor Green
