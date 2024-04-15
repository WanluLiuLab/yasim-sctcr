#!/usr/bin/env bash
set -ue
axel https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.refseq_rna.fna.gz
gunzip MANE.GRCh38.v1.3.refseq_rna.fna.gz
salmon index -t MANE.GRCh38.v1.3.refseq_rna.fna -i MANE.GRCh38.v1.3.refseq_rna.salmon_idx.d -p 40

# https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0
axel https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam
axel https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam.bai
axel https://cf.10xgenomics.com/samples/cell-vdj/7.0.1/SC5pv2_GEX_Human_Lung_Carcinoma_DTC/SC5pv2_GEX_Human_Lung_Carcinoma_DTC_possorted_genome_bam.bam
axel https://cf.10xgenomics.com/samples/cell-vdj/7.0.1/SC5pv2_GEX_Human_Lung_Carcinoma_DTC/SC5pv2_GEX_Human_Lung_Carcinoma_DTC_possorted_genome_bam.bam.bai

for sample_name in SC5pv2_GEX_Human_Lung_Carcinoma_DTC_possorted_genome_bam pbmc_1k_v3_possorted_genome_bam; do
    samtools flagstat -O json "${sample_name}".bam >"${sample_name}".bam.flagstat.json
    samtools fastq "${sample_name}".bam | pigz -9 -cf >"${sample_name}".fq.gz
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
done

# https://www.10xgenomics.com/datasets/nsclc-tumor-tcr-enrichment-from-amplified-c-dna-1-standard-2-2-0
axel https://cf.10xgenomics.com/samples/cell-vdj/2.2.0/vdj_v1_hs_nsclc_t/vdj_v1_hs_nsclc_t_all_contig.bam
axel https://cf.10xgenomics.com/samples/cell-vdj/2.2.0/vdj_v1_hs_nsclc_t/vdj_v1_hs_nsclc_t_all_contig.bam.bai
axel https://cf.10xgenomics.com/samples/cell-vdj/2.2.0/vdj_v1_hs_nsclc_t/vdj_v1_hs_nsclc_t_all_contig.fasta
axel https://cf.10xgenomics.com/samples/cell-vdj/2.2.0/vdj_v1_hs_nsclc_t/vdj_v1_hs_nsclc_t_all_contig_annotations.bed

axel https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_all_contig.bam
axel https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_all_contig.bam.bai
axel https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_all_contig.fasta
axel https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_all_contig_annotations.bed

for sample_name in sc5p_v2_hs_PBMC_10k_t_all_contig vdj_v1_hs_nsclc_t_all_contig; do
    samtools depth -@40 -aa -H "${sample_name}".bam |
        xz -9 -T0 -vvv >"${sample_name}".tcr.depth.tsv.xz
done
samtools fastq \
    -@ 40 \
    vdj_v1_hs_nsclc_t_all_contig.bam \
    -1 vdj_v1_hs_nsclc_t_all_contig_1.fq \
    -2 vdj_v1_hs_nsclc_t_all_contig_2.fq \
    -0 /dev/null \
    -s /dev/null \
    -N

fastp -i vdj_v1_hs_nsclc_t_all_contig_1.fq \
    -I vdj_v1_hs_nsclc_t_all_contig_2.fq \
    --stdout \
    --html vdj_v1_hs_nsclc_t_all_contig.fastp.html \
    --thread 40 >/dev/null
