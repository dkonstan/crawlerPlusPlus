#!/bin/bash

./crawler++ -t top.txt -c crd.xyz -i min.in -l log -o traj.xyz -x out.xyz;
./crawler++ -t top.txt -c out.xyz -i run.in -l log2 -o traj.xyz -x out.xyz;
