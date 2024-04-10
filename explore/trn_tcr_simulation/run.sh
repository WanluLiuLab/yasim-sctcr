#!/usr/bin/env bash
set -ue

python -m yasim art \
    -F sim_tcr.nt.fa.d \
    -o sim_tcr_50 \
    --sequencer_name GA2 \
    --read_length 50 \
    -d tcr_depth.tsv \
    -e art_illumina \
    -j 20
python -m yasim art \
    -F sim_tcr.nt.fa.d \
    -o sim_tcr_100 \
    --sequencer_name HS20 \
    --read_length 100 \
    -d tcr_depth.tsv \
    -e art_illumina \
    -j 20
python -m yasim art \
    -F sim_tcr.nt.fa.d \
    -o sim_tcr_150 \
    --sequencer_name HS25 \
    --read_length 150 \
    -d tcr_depth.tsv \
    -e art_illumina \
    -j 20
python -m yasim art \
    -F sim_tcr.nt.fa.d \
    -o sim_tcr_250 \
    --sequencer_name MSv3 \
    --read_length 250 \
    -d tcr_depth.tsv \
    -e art_illumina \
    -j 20
# See: sim_tcr.d

# Gene
python -m yasim_sc sample_pcg \
    -i data/ncbi_dataset.tsv \
    -g hg38.ncbiRefSeq.gtf \
    -o hg38.ncbiRefSeq_subsampled.gtf \
    --num_genes_to_sample 20
python -m yasim generate_gene_depth \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -d 5 \
    -o notcr_gene_depth.tsv
python -m yasim generate_isoform_depth \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -d notcr_gene_depth.tsv \
    -o notcr_isoform_depth.tsv
python -m yasim_sc generate_barcoded_isoform_replicates \
    -d notcr_isoform_depth.tsv \
    -b barcode.txt \
    -o notcr_isoform_depth_sc.d
rm -fr notcr_trans.fa.d
python -m labw_utils.bioutils transcribe \
    -f hg38.fa \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -o notcr_trans.fa

python -m yasim_sc art \
    -F notcr_trans.fa.d \
    -o notcr_trans_sim.d \
    -d notcr_isoform_depth_sc.d \
    -e art_illumina \
    -j 20

mkdir -p trust4_result

run-trust4 -u sim_tcr.fq \
    -f trust4_index/bcrtcr.fa \
    -t 40 \
    --ref trust4_index/IMGT+C.fa \
    --od trust4_result
