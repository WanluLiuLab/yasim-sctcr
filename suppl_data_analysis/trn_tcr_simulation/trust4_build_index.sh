#!/usr/bin/env bash
mkdir -p trust4_index
cd trust4_index || exit 1
axel https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz
axel https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
gunzip ./*.gz
perl "${CONDA_PREFIX}"/bin/BuildImgtAnnot.pl Homo_sapien >IMGT+C.fa
grep ">" IMGT+C.fa | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq >bcr_tcr_gene_name.txt

perl "${CONDA_PREFIX}"/bin/BuildDatabaseFa.pl \
    GRCh38.primary_assembly.genome.fa \
    gencode.v43.basic.annotation.gtf \
    bcr_tcr_gene_name.txt >bcrtcr.fa
