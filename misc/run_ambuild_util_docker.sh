#!/bin/bash

# Usage:
# ambuild_docker.py <volume_arguments> script

# Get root dir and script argumnts
ambuild_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"

# Where python lives in the container
docker_python_install_path=/usr/lib/python3/dist-packages

# Run
docker run \
-it \
--rm \
--runtime=nvidia \
--volume "$PWD":/home/glotzerlab \
--volume ${ambuild_dir}/ambuild:${docker_python_install_path}/ambuild \
--env PYTHONPATH=${docker_python_install_path} \
--workdir /home/glotzerlab \
glotzerlab/software \
python3 ${docker_python_install_path}/ambuild/ab_util.py $*

