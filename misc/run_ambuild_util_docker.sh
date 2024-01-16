#!/bin/bash

# Usage:
# ambuild_docker.py <volume_arguments> script

# Get root dir and script argumnts
ambuild_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"

# Run
docker run \
-it \
--rm \
--runtime=nvidia \
--volume "$PWD":/home/glotzerlab \
--volume ${ambuild_dir}/ambuild:/usr/lib/python3/dist-packages/ambuild \
--workdir /home/glotzerlab \
glotzerlab/software:2020.11.18-cuda10 \
python3 /usr/lib/python3/dist-packages/ambuild/ab_util.py $*



