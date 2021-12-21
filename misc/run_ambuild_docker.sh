#!/bin/bash

# Usage:
# ambuild_docker.py <volume_arguments> script

# Get root dir and script argumnts
run_dir="$PWD"
ambuild_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"
script=$1 
extra_args=""
if [ ${#} -ge 2 ]; then
   script="${@:(-1):1}"
   extra_args="${@:1:$(($#-1))}"
fi

# Run
docker run \
--rm \
--runtime=nvidia \
--volume $run_dir:$run_dir \
--volume ${ambuild_dir}/ambuild:/usr/lib/python3/dist-packages/ambuild \
--env PYTHONPATH=/usr/lib/python3/dist-packages \
--workdir $run_dir \
$extra_args \
glotzerlab/software \
python3 $script
