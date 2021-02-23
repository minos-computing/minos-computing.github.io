#!/bin/bash

num_instances=`ls nvdla_emulator_configs/ | wc -l`


devices=""


for i in `seq 0 $(( num_instances-1 ))`; do
devices="${devices}nvdlaemu "
port=$(( 6000+i )) 
export POCL_NVDLAEMU${i}_PARAMETERS="127.0.0.1:${port}"
done
POCL_DEVICES=${devices} mcl_sched &