#!/bin/bash

# This is somewhat of a hack to link in the test and blocks directories
# into the python directory so the test scripts can find them. It'll do
# for now but a better approach should be found as linking non-python
# directories into the python installation is a bad idea!

# Get root dir and script argumnts
run_dir="$PWD"
ambuild_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"

# Run
--runtime=nvidia \
docker run \
--rm \
--volume $run_dir:$run_dir \
--volume ${ambuild_dir}/ambuild:/usr/lib/python3/dist-packages/ambuild \
--volume ${ambuild_dir}/tests:/usr/lib/python3/dist-packages/tests \
--volume ${ambuild_dir}/blocks:/usr/lib/python3/dist-packages/blocks \
--workdir $run_dir \
glotzerlab/software \
python3 $*
