#!/bin/bash

run_dir=`pwd -P`
ambuild_dir=$(cd $run_dir/..;pwd -P)

# Get current user id
uid=`id | awk -F "=" {'print $2'}| awk -F "(" {'print $1'}`

args="python3 ./run_tests.py"
# Below for runnining individual unittests e.g: ./run_tests_docker.sh python3 -m unittest testCell.Test.testCat1Paf2
if [ $# -gt 0 ]; then
    args=$*
fi

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
--user $uid \
glotzerlab/software \
$args
