#!/usr/bin/env bash
set -ue

axel https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam
axel https://cf.10xgenomics.com/samples/cell-vdj/7.0.1/SC5pv2_GEX_Human_Lung_Carcinoma_DTC/SC5pv2_GEX_Human_Lung_Carcinoma_DTC_possorted_genome_bam.bam
axel https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.refseq_rna.fna.gz

salmon index -t MANE.GRCh38.v1.3.refseq_rna.fna.gz -i MANE.GRCh38.v1.3.refseq_rna.salmon_idx.d -p 40

sample_name=SC5pv2_GEX_Human_Lung_Carcinoma_DTC_possorted_genome_bam

samtools index "${sample_name}".bam
samtools flagstat "${sample_name}".bam >"${sample_name}".bam.flagstat.txt
samtools fastq "${sample_name}".bam | pigz -9 -T0 -cf >"${sample_name}".fq.gz
salmon --no-version-check quant \
    -i MANE.GRCh38.v1.3.refseq_rna.salmon_idx.d \
    -l U \
    -r "${sample_name}".fq.gz \
    -o "${sample_name}".salmon.d \
    --writeMappings \
    -p 40 | samtools view -o "${sample_name}".salmon.bam
samtools sort -@40 "${sample_name}".salmon.bam -o "${sample_name}".salmon.sorted.bam
samtools index "${sample_name}".salmon.sorted.bam
samtools depth -@40 -aa -H "${sample_name}".salmon.sorted.bam |
    xz -9 -T0 -vvv >"${sample_name}".salmon.sorted.depth.tsv.xz
