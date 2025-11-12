#!/bin/bash
#SBATCH -o snakePipe_%A.log
#SBATCH -e snakePipe_%A.err
#SBATCH -t 1-0
#SBATCH -J snakePipe
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=away@email.unc.edu

module load python

ROOT_DIR="$(pwd)"


export MAMBA_ROOT_PREFIX="$ROOT_DIR/.mamba_root"
export TMPDIR="$ROOT_DIR/tmp_snakemake"
export _JAVA_OPTIONS="-Djava.io.tmpdir=$TMPDIR"

rm -rf "$TMPDIR"
mkdir -p "$MAMBA_ROOT_PREFIX" "$TMPDIR"

snakemake \
  --directory "$ROOT_DIR" \
  --snakefile "$ROOT_DIR/Snakefile" \
  --printshellcmds \
  --jobs=100 \
  --latency-wait 300 \
  --rerun-incomplete \
  --cluster-config "$ROOT_DIR/cluster_config.json" \
  --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -n {cluster.ntasks} -o {cluster.output} -e {cluster.err} -J {cluster.jobname}"
