#!/usr/bin/env bash
#SBATCH --job-name=run_sim_tcr
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=48:00:00
# shellcheck disable=SC2002

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim_dev

cd /slurm/home/yrd/liulab/yuzhejian/yasim-sctcr
. bashrc
cd /slurm/home/yrd/liulab/yuzhejian/yasim-sctcr/explore/trn_tcr_simulation
set -ue

data_name=HU_0043_Blood_10x

for tcr_depth in 2 4 6 8 10 20 40 60 80 100; do
    python -m yasim art \
        -F sim/"${data_name}".sim.d/sim_t_cell.rc.nt.fa.d \
        -o sim/"${data_name}".sim.d/art_diff_depth/sim_tcr_rlen250_tcrd${tcr_depth} \
        --sequencer_name MSv3 \
        --read_length 250 \
        -d sim/"${data_name}".sim.d/scTCR.depth"${tcr_depth}".tsv \
        -e art_illumina \
        -j 40 \
        --amplicon
    rm -f sim/"${data_name}".sim.d/art_diff_depth/sim_tcr_rlen250_tcrd${tcr_depth}.fq
done

RLENS=(50 100 150 250)
MODELS=(GA2 HS20 HS25 MSv3)
tcr_depth=400
for i in {0..3}; do
    python -m yasim art \
        -F sim/"${data_name}".sim.d/sim_t_cell.rc.nt.fa.d \
        -o sim/"${data_name}".sim.d/art_diff_rlen/sim_tcr_rlen"${RLENS[i]}"_tcrd400 \
        --sequencer_name "${MODELS[i]}" \
        --read_length "${RLENS[i]}" \
        -d sim/"${data_name}".sim.d/scTCR.depth"${tcr_depth}".tsv \
        -e art_illumina \
        -j 40 \
        --amplicon
    rm -f sim/"${data_name}".sim.d/art_diff_rlen/sim_tcr_rlen"${RLENS[i]}"_tcrd400.fq
done

for t_cell_num in 100 500 1000; do
    for replicate_num in {1..10}; do
        prefix="${data_name}"_sim_tcell_only_ncells"${t_cell_num}"_rep"${replicate_num}"
        python -m yasim art \
            -F sim/"${prefix}".d/sim_t_cell.rc.nt.fa.d \
            -o sim/"${prefix}".d/art_sim_t_cell_rlen250 \
            --sequencer_name MSv3 \
            --read_length 250 \
            -d sim/"${prefix}".d/scTCR.depth10.tsv \
            -e art_illumina \
            -j 40 \
            --amplicon
        rm -f sim/"${prefix}".d/art_sim_t_cell_rlen250.fq
    done
done
