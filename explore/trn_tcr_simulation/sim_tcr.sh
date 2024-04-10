#!/usr/bin/env bash
python -m \

python -m yasim_sctcr rearrange_tcr \
    --tcr_cache_path ref/tcr_cache.json \
    --cdr3_deletion_table_path data/cdr3_deletion_table.json \
    --cdr3_insertion_table_path data/cdr3_insertion_table.json \
    --usage_bias_json data/usage_bias.json \
    -n "$(cat sim/HU_0043_Blood_10x.sim.d/t_cell_bc.txt | wc -l)" \
    -o sim/HU_0043_Blood_10x.sim.d/sim_tcr
python -m yasim_sctcr generate_tcr_clonal_expansion \
    -b sim/HU_0043_Blood_10x.sim.d/t_cell_bc.txt \
    --src_tcr_stats_tsv sim/HU_0043_Blood_10x.sim.d/sim_tcr.stats.tsv \
    --dst_nt_fasta sim/HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa \
    --alpha 1
python -m labw_utils.bioutils split_fasta sim/HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa

for tcr_depth in 2 4 8 10 20 40 60 80 100; do
    python -m yasim_sctcr generate_tcr_depth \
        -b sim/HU_0043_Blood_10x.sim.d/t_cell_bc.txt \
        -o sim/HU_0043_Blood_10x.sim.d/tcr_depth.d/"${tcr_depth}".tsv \
        -d "${tcr_depth}"
    python -m yasim art \
        -F sim/HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa.d \
        -o sim/HU_0043_Blood_10x.sim.d/diff_depth/sim_tcr_rlen150_tcrd${tcr_depth} \
        --sequencer_name HS25 \
        --read_length 150 \
        -d sim/HU_0043_Blood_10x.sim.d/tcr_depth.d/"${tcr_depth}".tsv \
        -e art_illumina \
        -j 40
    rm -f sim/HU_0043_Blood_10x.sim.d/diff_depth/sim_tcr_rlen150_tcrd${tcr_depth}.fq
done

RLENS=(50 100 150 250)
MODELS=(GA2 HS20 HS25 MSv3)
for i in {0..3}; do
    python -m yasim art \
        -F sim/HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa.d \
        -o sim/HU_0043_Blood_10x.sim.d/diff_rlen/sim_tcr_rlen"${RLENS[i]}"_tcrd400 \
        --sequencer_name "${MODELS[i]}" \
        --read_length "${RLENS[i]}" \
        -d sim/HU_0043_Blood_10x.sim.d/tcr_depth.d/400.tsv \
        -e art_illumina \
        -j 40
    rm -f sim/HU_0043_Blood_10x.sim.d/diff_rlen/sim_tcr_rlen"${RLENS[i]}"_tcrd400.fq
done

NTCELLS=(100 500 1000)



# 1. No TCR clonal expansion, scRNA (1.0)+scTCR (depth=2, 4, 6, 8, 10, 20, 40, 60, 80, 100, len=150SE) + readLen (50, 100, 150, 250SE, TCR depth=400).
# 2. (T-cell only) 100, 500, 1000 cells * 10 sample. Out: Clonal expansion in each sample. depth: scRNA: 0.8, scTCR: 10X
