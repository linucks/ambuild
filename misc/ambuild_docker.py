#!/bin/bash
/usr/local/bin/docker run \
-it \
--rm \
--volume "$PWD":/home/glotzerlab \
--volume /opt/ambuild/ambuild:/usr/lib/python3/dist-packages/ambuild \
--volume /opt/ambuild/params:/home/glotzerlab/params \
--volume /opt/ambuild/blocks:/home/glotzerlab/blocks \
--volume /opt/ambuild/tests/test_data:/home/glotzerlab/test_data \
--workdir /home/glotzerlab \
glotzerlab/software \
python3 $*

#--volume $HOME:$HOME \
#--runtime=nvidia \
#--volume $HOME/Dropbox/Ambuild_Files/Parameters:$HOME/Dropbox/Ambuild_Files/Parameters \
#--volume $HOME/Dropbox/Ambuild_Files/Blocks:$HOME/Dropbox/Ambuild_Files/Blocks \
