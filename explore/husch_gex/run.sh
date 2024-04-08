#!/usr/bin/env bash
set -ue
mkdir -p raw_data
cd raw_data
# 10X: HU_0043_Blood_10x
# Indrop: HU_0196_Kidney_GSE109564
# Smart-Seq2: HU_0148_Decidua_EBI
# Microwell: HU_0223_Muscle_GSE134355
# Drop-Seq: HU_0125_Cerebrospinal-Fluid_GSE134577

mkdir -p figs parquets

for data_name in \
    HU_0196_Kidney_GSE109564 \
    HU_0043_Blood_10x \
    HU_0223_Muscle_GSE134355 \
    HU_0148_Decidua_EBI \
    HU_0125_Cerebrospinal-Fluid_GSE134577; do
    if [ ! -f "${data_name}"_gene_count.h5 ]; then
        aws s3 cp --no-sign-request \
            s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_gene_count.h5 .
        aws s3 cp --no-sign-request \
            s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_meta.txt .
    fi
done
data_name=HU_0142_Breast_GSE138536
aws s3 cp --no-sign-request s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_gene_count.h5 sample_gene_count.h5
aws s3 cp --no-sign-request s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_meta.txt sample_meta.txt

mkdir -p ref
cd ref
axel https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2024-03-01.tsv
axel https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.ensembl_rna.fna.gz
gunzip MANE.GRCh38.v1.3.ensembl_rna.fna.gz
cat hgnc_complete_set_2024-03-01.tsv |
    grep -v '\t$' |
    cut -f 2,53 |
    sed 's;";;' |
    cut -f 1 -d '|' |
    sed -E 's;(\S+)\s(\S+);\2\t\1;' |
    sed '1d' >hgnc_salmon_genemap.tsv
cd ..

for data_name in \
    HU_0196_Kidney_GSE109564 \
    HU_0043_Blood_10x \
    HU_0223_Muscle_GSE134355 \
    HU_0148_Decidua_EBI \
    HU_0125_Cerebrospinal-Fluid_GSE134577; do
    python -m yasim_sctcr scaffold \
        --transcript_gene_mapping ref/hgnc_salmon_genemap.tsv \
        --src_sc_data parquets/"${data_name}"_sim.parquet \
        --out "${data_name}".sim.d \
        --t_cell_regex 'CD[48]T'
done
