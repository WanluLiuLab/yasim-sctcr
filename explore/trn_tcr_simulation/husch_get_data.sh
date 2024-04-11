#!/usr/bin/env bash
# shellcheck disable=SC2002
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
    HU_0043_Blood_10x \
    HU_0196_Kidney_GSE109564 \
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

Rscript husch_qc.R
