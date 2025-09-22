param()

# Generate pinned requirements.txt from requirements.in using pip-compile
# Requires pip-tools: pip install pip-tools

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
$projectRoot = Join-Path $scriptDir '..' | Resolve-Path
Set-Location $projectRoot

if (-not (Get-Command pip-compile -ErrorAction SilentlyContinue)) {
    Write-Host "pip-compile not found. Installing pip-tools..."
    pip install --user pip-tools
}

Write-Host "Compiling requirements.txt from requirements.in..."
pip-compile --output-file=requirements.txt requirements.in
Write-Host "Wrote requirements.txt"
