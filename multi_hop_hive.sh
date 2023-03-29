#!/bin/bash
if [ $# -lt 1 ]; then
  echo "usage: multihop host port"
else
  host1=jbrandt35@login-hive.pace.gatech.edu
  host2=$1
  port=$2
  ssh -L ${port}:localhost:${port} $host1 \
    ssh -L ${port}:localhost:${port} -N $host2
fi

