#!/bin/bash
N_Cores=$1
docker run --privileged --rm \
--mount type=bind,source="$(pwd)"/../config,target=/config \
--mount type=bind,source="$(pwd)"/benchmarks,target=/home/benchmarks \
--mount type=bind,source="$(pwd)"/results,target=/home/results \
--mount type=bind,source="$(pwd)"/../raw_data,target=/raw_data \
--mount type=bind,source="$(pwd)"/scripts,target=/home/scripts zmxu/g2g_snakemake_docker snakemake --cores $N_Cores --rerun-incomplete 
