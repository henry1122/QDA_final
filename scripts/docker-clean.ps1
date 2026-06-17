# Free Docker disk space after heavy qsyn builds.
# Restart Docker Desktop first if you see API 500 errors.
Write-Host "Docker disk usage (before):" -ForegroundColor Cyan
docker system df 2>$null

Write-Host "`nRemoving stopped containers..." -ForegroundColor Yellow
docker container prune -f

Write-Host "Removing dangling images..." -ForegroundColor Yellow
docker image prune -f

Write-Host "Removing build cache (this frees the most space)..." -ForegroundColor Yellow
docker builder prune -a -f

Write-Host "`nOptional: remove qsyn test images (uncomment if needed)"
Write-Host "  docker rmi qsyn-test-gcc qsyn-local 2>`$null"

Write-Host "`nDocker disk usage (after):" -ForegroundColor Cyan
docker system df
