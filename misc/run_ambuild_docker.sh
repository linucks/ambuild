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
--volume ${ambuild_dir}/ambuild:/home/abbie/ambuild/ambuild \
--workdir $run_dir \
$extra_args \
glotzerlab/software:2020.11.18-cuda10 \
python3 $script
