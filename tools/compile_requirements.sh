This file has been removed as the Windows batch wrapper exists and has been tested.
#!/usr/bin/env bash
set -euo pipefail

# Generate pinned requirements.txt from requirements.in using pip-compile
# Requires: pip install pip-tools

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"

if ! command -v pip-compile >/dev/null 2>&1; then
  echo "pip-compile not found. Installing pip-tools..."
  pip install --user pip-tools
fi

echo "Compiling requirements.txt from requirements.in..."
pip-compile --output-file=requirements.txt requirements.in
echo "Wrote requirements.txt"
