# Run pp+ vs TODD+PP benchmark in Docker. Writes results/benchmark_*_pp.txt and *_todd_pp.txt
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)

$image = "qsyn-test-gcc"
if (-not (docker image inspect $image 2>$null)) {
    Write-Host "Building $image..." -ForegroundColor Cyan
    docker build -t $image -f docker/gcc-test.Dockerfile .
}

Write-Host "Building benchmark-phasepoly + running pp+ vs TODD+PP..." -ForegroundColor Cyan
docker run --rm `
    --entrypoint /bin/bash `
    --security-opt seccomp=unconfined `
    -v "${PWD}:/app/qsyn" `
    -v qsyn-build-cache:/tmp/qsyn-build `
    $image `
    -lc "cp /app/qsyn/scripts/docker-todd-pp-bench.sh /tmp/run.sh && sed -i 's/\r$//' /tmp/run.sh && bash /tmp/run.sh"

if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Done. Check results/benchmark_*_pp.txt and results/benchmark_*_todd_pp.txt" -ForegroundColor Green
