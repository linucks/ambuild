#!/bin/bash

[[ -f trajectory.xyz ]] && rm trajectory.xyz
cat `ls *.xyz | grep -v "_P" | sort -t "_" -k 2.1  -n ` > trajectory.xyz
