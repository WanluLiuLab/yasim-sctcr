#!/usr/bin/env bash
set -ue
mkdir -p trust4_index
cd trust4_index || exit 1
axel -4 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz
axel -4 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
gunzip ./*.gz
perl "${CONDA_PREFIX}"/bin/BuildImgtAnnot.pl Homo_sapien >IMGT+C.fa
grep ">" IMGT+C.fa | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq >bcr_tcr_gene_name.txt

perl "${CONDA_PREFIX}"/bin/BuildDatabaseFa.pl \
    GRCh38.primary_assembly.genome.fa \
    gencode.v43.basic.annotation.gtf \
    bcr_tcr_gene_name.txt >bcrtcr.fa
cd ..

python -m rearrange_tcr

python -m yasim_sctcr rearrange_tcr \
    --tcr_cache_path /home/yuzj/Documents/yasim-sctcr/explore/trn_tcr_simulation/ref/tcr_cache.json \
    --cdr3_deletion_table_path ../data/cdr3_deletion_table.min.json.xz \
    --cdr3_insertion_table_path ../data/cdr3_insertion_table.min.json.xz \
    --usage_bias_json ../data/usage_bias.min.json.xz \
    -n 100 \
    -o sim_tcr
python -m yasim_sc _generate_barcode \
    -n 100 \
    -l 15 \
    -o bc.txt
python -m yasim_sctcr generate_tcr_clonal_expansion \
    -b bc.txt \
    --src_tcr_stats_tsv sim_tcr.stats.tsv \
    --dst_nt_fasta sim_t_cell.nt.fa \
    --alpha 1
python -m yasim_sctcr generate_tcr_depth \
    -b bc.txt \
    -o scTCR.depth.tsv \
    -d 20
seqkit seq --reverse <sim_t_cell.nt.fa >sim_t_cell.rev.nt.fa
seqkit seq --reverse --complement <sim_t_cell.nt.fa >sim_t_cell.rc.nt.fa

python -m labw_utils.bioutils split_fasta sim_t_cell.nt.fa
python -m labw_utils.bioutils split_fasta sim_t_cell.rc.nt.fa
python -m labw_utils.bioutils split_fasta sim_t_cell.rev.nt.fa

python -m yasim art \
    -F sim_t_cell.nt.fa.d \
    -o sim_tcr_batch \
    --sequencer_name HS25 \
    --read_length 150 \
    -d scTCR.depth.tsv \
    -e art_illumina \
    --preserve_intermediate_files \
    -j 20

python -m yasim art \
    -F sim_t_cell.nt.fa.d \
    -o sim_tcr_ampli \
    --sequencer_name HS25 \
    --read_length 150 \
    -d scTCR.depth.tsv \
    -e art_illumina \
    --preserve_intermediate_files \
    -j 20 \
    --amplicon

python -m yasim art \
    -F sim_t_cell.rc.nt.fa.d \
    -o sim_tcr_ampli_rc \
    --sequencer_name HS25 \
    --read_length 150 \
    -d scTCR.depth.tsv \
    -e art_illumina \
    --preserve_intermediate_files \
    -j 20 \
    --amplicon

python -m yasim art \
    -F sim_t_cell.rev.nt.fa.d \
    -o sim_tcr_ampli_rev \
    --sequencer_name HS25 \
    --read_length 150 \
    -d scTCR.depth.tsv \
    -e art_illumina \
    --preserve_intermediate_files \
    -j 20 \
    --amplicon

mkdir -p trust4_result
for name in sim_tcr_batch sim_tcr_ampli sim_tcr_ampli_rc sim_tcr_ampli_rev; do
    run-trust4 \
        -u "${name}".fq \
        -t 40 \
        -f trust4_index/IMGT+C.fa \
        --od trust4_result/"${name}" &>trust4_result/"${name}".log
    wc -l trust4_result/"${name}"/TRUST_"${name}"_report.tsv
done
# /home/yuzj/Documents/yasim-sctcr/test/trust4_result/

run-trust4 \
    -u sim_tcr.fq \
    -f trust4_index/bcrtcr.fa \
    -t 40 \
    --ref trust4_index/IMGT+C.fa \
    --od trust4_result
