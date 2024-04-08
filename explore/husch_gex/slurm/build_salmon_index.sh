#!/usr/bin/env bash
#SBATCH --job-name=build_salmon_index
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks-per-node=4
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=48:00:00


eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim-salmon
set -ue

cd /slurm/home/yrd/liulab/yuzhejian/yasim-sctcr/explore/husch_gex

salmon index -t MANE.GRCh38.v1.3.refseq_rna.fna -i MANE.GRCh38.v1.3.refseq_rna.salmon_idx -p 4
