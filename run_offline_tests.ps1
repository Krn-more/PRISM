$ErrorActionPreference = 'Stop'
$runtime = Join-Path $PSScriptRoot '.build-venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $runtime)) {
    throw 'Bundled build virtual environment is unavailable. Create it before running offline tests.'
}

# These tests are unit/regression tests only. They must make no network calls.
& $runtime -m unittest discover -p 'test_*.py'
exit $LASTEXITCODE
