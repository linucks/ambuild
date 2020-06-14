#!/bin/bash
docker run \
-it \
--rm \
--runtime=nvidia \
--volume "$PWD":/home/glotzerlab \
--volume /opt/ambuild/ambuild:/usr/lib/python3/dist-packages/ambuild \
--volume /opt/ambuild/tests/params:/home/glotzerlab/params \
--volume /opt/ambuild/tests/blocks:/home/glotzerlab/blocks \
--workdir /home/glotzerlab \
glotzerlab/software \
python3 $*

#--volume $HOME:$HOME \
#--runtime=nvidia \
#--volume $HOME/Dropbox/Ambuild_Files/Parameters:$HOME/Dropbox/Ambuild_Files/Parameters \
#--volume $HOME/Dropbox/Ambuild_Files/Blocks:$HOME/Dropbox/Ambuild_Files/Blocks \
