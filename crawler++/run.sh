#!/bin/bash

./crawler++ -t top.txt -c crd.xyz -i min.in -l log -o traj1.xyz -x out.xyz;
./crawler++ -t top.txt -c out.xyz -i run.in -l log2 -o traj2.xyz -x out2.xyz;


./crawler++ -t Li+W20_tip3p.top -c Li_waterClusterReordered.xyz -i min.in -l log -o traj1.xyz -x out.xyz;

./crawler++ -t Li+W20_tip3p.top -c out.xyz -i run.in -l log -o traj2.xyz -x out2.xyz;