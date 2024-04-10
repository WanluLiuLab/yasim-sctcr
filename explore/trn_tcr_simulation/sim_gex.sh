#!/usr/bin/env bash
set -ue
Rscript husch_sim.R

data_name=HU_0043_Blood_10x

python -m yasim_sctcr scaffold \
    --transcript_gene_mapping ref/ens.trans_gene_map.tsv \
    --src_sc_data parquets/"${data_name}"_sim.parquet \
    --out sim/"${data_name}".sim.d \
    --t_cell_regex 'CD[48]T'

for t_cell_num in 100 500 1000; do
    python -m yasim_sctcr scaffold \
        --transcript_gene_mapping ref/ens.trans_gene_map.tsv \
        --src_sc_data parquets/"${data_name}"_sim_tcell_only_"${t_cell_num}".parquet \
        --out sim/"${data_name}".sim.tcell_only_"${t_cell_num}".d \
        --t_cell_regex 'CD[48]T'
done
