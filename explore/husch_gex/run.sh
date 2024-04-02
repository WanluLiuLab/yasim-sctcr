#!/usr/bin/env bash
set -ue
mkdir -p raw_data
cd raw_data
for data_name in HU_0196_Kidney_GSE109564 HU_0281_Blood_GSE138867 HU_0133_Cerebral-Cortex_GSE134355; do
    aws s3 cp --no-sign-request s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_gene_count.h5 .
    aws s3 cp --no-sign-request s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_meta.txt .
done
data_name=HU_0142_Breast_GSE138536
aws s3 cp --no-sign-request s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_gene_count.h5 sample_gene_count.h5
aws s3 cp --no-sign-request s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_meta.txt sample_meta.txt
axel https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

mkdir -p gex_sim
mv ./*.csv gex_sim
tar cvf gex_sim.tar gex_sim
xz -9 -T0 -vv gex_sim.tar
