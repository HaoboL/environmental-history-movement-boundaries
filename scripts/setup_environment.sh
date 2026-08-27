#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
venv_dir="$repo_root/.venv"
python_command="${PYTHON:-python3}"

echo "[1/5] Checking Python"
"$python_command" -c 'import sys; assert sys.version_info >= (3, 10), "Python 3.10 or newer is required"; print(sys.version)'

echo "[2/5] Creating isolated environment at $venv_dir"
"$python_command" -m venv "$venv_dir"

echo "[3/5] Installing analysis and download dependencies"
"$venv_dir/bin/python" -m pip install --upgrade pip setuptools wheel
"$venv_dir/bin/python" -m pip install -r "$repo_root/requirements.txt"

echo "[4/5] Installing the radial-drawdown movement-segmentation dependency"
"$venv_dir/bin/python" -m pip install "git+https://github.com/HaoboL/radial-drawdown.git"

echo "[5/5] Running release validation and unit tests"
cd "$repo_root"
"$venv_dir/bin/python" scripts/validate_release.py
"$venv_dir/bin/python" -m unittest discover -s tests -v

echo "SETUP_COMPLETE"
echo "Activate later with: source .venv/bin/activate"
