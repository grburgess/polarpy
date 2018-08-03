#!/usr/bin/env bash
set -e

cd root/polarpy

pip install pytest 

python -m pytest -vv --cov=polarpy
