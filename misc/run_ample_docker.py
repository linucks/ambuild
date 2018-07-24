#!/bin/bash
docker run --runtime=nvidia --rm \
--volume "$PWD":/home/glotzerlab \
--volume /opt/ambuild_v1.0:/opt/ambuild_v1.0 \
--volume /home/pierre/Dropbox/Ambuild_Files/Parameters:/home/pierre/Dropbox/Ambuild_Files/Parameters \
--volume /home/pierre/Dropbox/Ambuild_Files/Blocks:/home/pierre/Dropbox/Ambuild_Files/Blocks \
-w /home/glotzerlab glotzerlab/software \
python3 /home/glotzerlab/$1
