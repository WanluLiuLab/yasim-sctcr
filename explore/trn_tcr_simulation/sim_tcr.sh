#!/usr/bin/env bash
set -ue
# shellcheck disable=SC2002

data_name=HU_0043_Blood_10x
python -m yasim_sctcr rearrange_tcr \
    --tcr_cache_path ref/tcr_cache.json \
    --cdr3_deletion_table_path data/cdr3_deletion_table.json \
    --cdr3_insertion_table_path data/cdr3_insertion_table.json \
    --usage_bias_json data/usage_bias.json \
    -n "$(cat sim/"${data_name}".sim.d/t_cell_bc.txt | wc -l)" \
    -o sim/"${data_name}".sim.d/sim_tcr
python -m yasim_sctcr generate_tcr_clonal_expansion \
    -b sim/"${data_name}".sim.d/t_cell_bc.txt \
    --src_tcr_stats_tsv sim/"${data_name}".sim.d/sim_tcr.stats.tsv \
    --dst_nt_fasta sim/"${data_name}".sim.d/sim_t_cell.nt.fa \
    --alpha 1
python -m labw_utils.bioutils split_fasta sim/"${data_name}".sim.d/sim_t_cell.nt.fa


for tcr_depth in 2 4 8 10 20 40 60 80 100; do
    python -m yasim_sctcr generate_tcr_depth \
        -b sim/"${data_name}".sim.d/t_cell_bc.txt \
        -o sim/"${data_name}".sim.d/scTCR.depth"${tcr_depth}".tsv \
        -d "${tcr_depth}"
done

NTCELLS=(100 500 1000)
NTCRS=(10 50 100)

for i in {0..2}; do
    t_cell_num=${NTCELLS[i]}
    for replicate_num in {1..10}; do
        prefix="${data_name}"_sim_tcell_only_ncells"${t_cell_num}"_rep"${replicate_num}"
        python -m yasim_sctcr rearrange_tcr \
            --tcr_cache_path ref/tcr_cache.json \
            --cdr3_deletion_table_path data/cdr3_deletion_table.json \
            --cdr3_insertion_table_path data/cdr3_insertion_table.json \
            --usage_bias_json data/usage_bias.json \
            -n "${NTCRS[i]}" \
            -o sim/"${prefix}".d/sim_tcr
        python -m yasim_sctcr generate_tcr_clonal_expansion \
            -b sim/"${prefix}".d/t_cell_bc.txt \
            --src_tcr_stats_tsv sim/"${prefix}".d/sim_tcr.stats.tsv \
            --dst_nt_fasta sim/"${prefix}".d/sim_t_cell.nt.fa \
            --alpha 1
        python -m labw_utils.bioutils split_fasta sim/"${prefix}".d/sim_t_cell.nt.fa
        python -m yasim_sctcr generate_tcr_depth \
                -b sim/"${prefix}".d/t_cell_bc.txt \
                -o sim/"${prefix}".d/scTCR.depth10.tsv \
                -d 10
    done
done

# 1. No TCR clonal expansion, scRNA (1.0)+scTCR (depth=2, 4, 6, 8, 10, 20, 40, 60, 80, 100, len=150SE) + readLen (50, 100, 150, 250SE, TCR depth=400).
# 2. (T-cell only) 100, 500, 1000 cells * 10 sample. Out: Clonal expansion in each sample. depth: scRNA: 0.8, scTCR: 10X