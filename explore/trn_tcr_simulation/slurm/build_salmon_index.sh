#!/usr/bin/env bash
#SBATCH --job-name=build_salmon_index
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks-per-node=10
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=48:00:00

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim-salmon
set -ue

cd /slurm/home/yrd/liulab/yuzhejian/yasim-sctcr/explore/trn_tcr_simulation/ref

salmon index -t ens.sel_genes.cdna.fa -i ens.sel_genes.cdna.fa.salmon_idx.d -p 10
