#!/usr/bin/env bash
set -e


# Make matplotlib non-interactive (otherwise it will crash
# all the tests)
export MPLBACKEND='Agg'
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd ~
rm -rf my_work_dir
mkdir my_work_dir
# Copy also dot files (.*)
shopt -s dotglob
cp -R ${TRAVIS_BUILD_DIR}/* my_work_dir/

cd my_work_dir

pip install pytest coverage pytest-cov codecov

python setup.py install

python -m pytest -vv --cov=polarpy

codecov -t 4ffa7470-e124-4021-a40c-6a6be5a560a4
