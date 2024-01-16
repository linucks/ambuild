  GNU nano 4.8                                                                           run_ambuild_docker_restart.sh                                                                                     
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

MAXLOOP=100

LOOPCOUNT=0
RESTART=1
while [ $RESTART -eq 1 ] && [ $LOOPCOUNT -lt $MAXLOOP ]; do
  echo $MAXLOOP
  echo Restart file is:
  ls -1t step_*.pkl.gz|head -1
  RFILE=`ls -1t step_*.pkl.gz|head -1`
  #STDIR=${> restart.o.$LOOPCOUNT 2> restart.e.$LOOPCOUNT}
  STDIR=" > restart.o.${LOOPCOUNT} 2> restart.e.${LOOPCOUNT}"
  FULLSTRING="${script} -i ${RFILE}"
  # Run
  docker run \
  --rm \
  --runtime=nvidia \
  --volume $run_dir:$run_dir \
  --volume ${ambuild_dir}/ambuild:/usr/lib/python3/dist-packages/ambuild \
  --workdir $run_dir \
  $extra_args \
  glotzerlab/software:2020.11.18-cuda10 \
  python3 $FULLSTRING > restart.o.${LOOPCOUNT} 2> restart.e.${LOOPCOUNT}
  mv runmd.log runmd.log.$LOOPCOUNT
  ./testrestart restart.e.$LOOPCOUNT
  RESTART=$?
  LOOPCOUNT=$[$LOOPCOUNT+1]
done

if [ $LOOPCOUNT -eq $MAXLOOP ]; then
   echo Error: Maximum restart loop count value of $MAXLOOP reached
fi


