#!/bin/bash
N_Cores=$0
snakemake --use-singularity --cores $N_Cores --rerun-incomplete --singularity-args "-B ./results:/home/results,../raw_data:/raw_data,./scripts:/home/scripts"
