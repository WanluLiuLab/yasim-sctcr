#!/usr/bin/env bash
set -ue
mkdir -p ref raw_data figs parquets
cd ref
axel https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
axel https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.pep.all.fa.gz
for name in cdna pep; do
    seqkit grep \
        --by-name \
        --use-regexp \
        -p 'chromosome:GRCh38:(7|14):' \
        Homo_sapiens.GRCh38."${name}".all.fa |
        seqkit grep \
            --by-name \
            --use-regexp \
            -p 'gene_biotype:TR_[VDJC]_gene' \
            /dev/stdin |
        seqkit grep \
            --by-name \
            --use-regexp \
            -p 'gene_symbol:TR[AB]' \
            /dev/stdin | sed -E 's;^>.+ gene_symbol:(\S+) .+$;>\1;' \
        >ens."${name}".fa
    samtools faidx ens."${name}".fa
done
seqtk subseq Homo_sapiens.GRCh38.cdna.all.fa <(cut -f 1 ens.trans_gene_map.tsv) >ens.sel_genes.cdna.fa

cd ..

python -m yasim_sctcr generate_tcr_cache \
    --tcr_cdna_fa_path ref/ens.cdna.fa \
    --tcr_pep_fa_path ref/ens.pep.fa \
    -o ref/tcr_cache.json

cd raw_data

data_name=HU_0043_Blood_10x
aws s3 cp --no-sign-request \
    s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_gene_count.h5 .
aws s3 cp --no-sign-request \
    s3://biostorage/HUSCH/HUSCH_data/"${data_name}"/"${data_name}"_meta.txt .
cd ..
Rscript run.R
