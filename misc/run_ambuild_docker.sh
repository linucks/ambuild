#!/bin/bash

# Usage:
# ambuild_docker.py <volume_arguments> script

# Get root dir and script argumnts
ambuild_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"
script=$1 
extra_args=""
if [ ${#} -ge 2 ]; then
   script="${@:(-1):1}"
   extra_args="${@:1:$(($#-1))}"
fi

# Run
docker run \
-it \
--rm \
--runtime=nvidia \
--volume "$PWD":/home/glotzerlab \
--volume ${ambuild_dir}/ambuild:/usr/lib/python3/dist-packages/ambuild \
--workdir /home/glotzerlab \
$extra_args \
glotzerlab/software \
python3 $script

