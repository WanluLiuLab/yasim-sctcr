#!/usr/bin/env bash
set -ue
mkdir -p raw_data
cd raw_data
# 10X: HU_0043_Blood_10x
# Indrop: HU_0196_Kidney_GSE109564
# Smart-Seq2: HU_0148_Decidua_EBI
# Microwell: HU_0223_Muscle_GSE134355
# Drop-Seq: HU_0125_Cerebrospinal-Fluid_GSE134577

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
axel https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

# https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set-03-01.txt

curl https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.refseq_rna.fna.gz | gunzip -cdf > MANE.GRCh38.v1.3.refseq_rna.fna
# curl https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.refseq_genomic.gtf.gz | gunzip -cdf > MANE.GRCh38.v1.3.refseq_genomic.gtf
curl https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.summary.txt.gz | gunzip > MANE.GRCh38.v1.3.summary.txt
cat MANE.GRCh38.v1.3.summary.txt | cut -f 4,6 | sed -E 's;(\S+)\s(\S+);\2\t\1;' | sed '1d' > MANE_salmon_genemap.tsv

