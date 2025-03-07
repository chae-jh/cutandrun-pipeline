#!/bin/bash

#SBATCH -o snakePipe_%A.log
#SBATCH -e snakePipe_%A.err
#SBATCH -t 1-0
#SBATCH -J snakePipe
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=away@email.unc.edu

snakemake --printshellcmds --jobs=100 --latency-wait 300 --cluster-config cluster_config.json --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -n {cluster.ntasks} -o {cluster.output} -e {cluster.err} -J {cluster.jobname}"
