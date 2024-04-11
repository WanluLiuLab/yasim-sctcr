#!/usr/bin/env bash
#SBATCH --job-name=run_simulation
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=40GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=48:00:00

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim_dev
set -ue
cd /slurm/home/yrd/liulab/yuzhejian/yasim-sctcr
. bashrc

cd /slurm/home/yrd/liulab/yuzhejian/yasim-sctcr/explore/trn_tcr_simulation

python scART.py
