#!/usr/bin/env bash
set -e

cd ~

pip install pytest 

python -m pytest -vv --cov=polarpy
