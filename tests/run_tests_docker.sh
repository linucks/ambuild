#!/bin/bash

run_dir=`pwd -P`
ambuild_dir=$(cd $run_dir/..;pwd -P)

# Not sure if setting PYTHONPATH below is sensible. In the container it's set as
# "/ignore/pythonpath" so I'm not sure if this is used by HOOMDBLUE internally
# However it seems to work and deals with the problem that different container versions
# use different python versions, which have their python packages stored in different locations

# Drop this for the time being as we won't always have this present
#--runtime=nvidia \

docker run \
--rm \
--volume $run_dir:$run_dir \
--volume ${ambuild_dir}/ambuild:/usr/lib/python3/dist-packages/ambuild \
--workdir $run_dir \
--env PYTHONPATH=/usr/lib/python3/dist-packages \
-it ubuntu /bin/bash ./script.sh

#glotzerlab/software \
#python3 ./run_tests.py
