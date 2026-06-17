# Run Table 1 benchmark suite in Docker (reuses qsyn-build-cache).
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)

$image = "qsyn-test-gcc"
if (-not (docker image inspect $image 2>$null)) {
    Write-Host "Building $image..." -ForegroundColor Cyan
    docker build -t $image -f docker/gcc-test.Dockerfile .
}

Write-Host "Building qsyn (if needed) + Table 1 benchmark..." -ForegroundColor Cyan
docker run --rm `
    --entrypoint /bin/bash `
    --security-opt seccomp=unconfined `
    -v "${PWD}:/app/qsyn:ro" `
    -v qsyn-build-cache:/tmp/qsyn-build `
    $image `
    -lc "cp /app/qsyn/scripts/table1-bench.sh /tmp/table1-bench.sh && sed -i 's/\r$//' /tmp/table1-bench.sh /app/qsyn/docker/test-entrypoint.sh 2>/dev/null; bash /app/qsyn/docker/test-entrypoint.sh '[phasepoly]' >/dev/null 2>&1; bash /tmp/table1-bench.sh"

if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
