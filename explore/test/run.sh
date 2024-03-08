#!/usr/bin/env bash
set -ue
axel https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory.zip
unzip IMGT_V-QUEST_reference_directory.zip
axel https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP
axel https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP
axel https://www.imgt.org/download/GENE-DB/IMGTGENEDB-GeneList
axel https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
axel https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
axel https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
mkdir blastdb
cd blastdb
wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar
tar xvf ncbi_human_c_genes.tar
wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/airr/airr_c_human.tar
tar xvf airr_c_human.tar
cd ..
wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz
tar xvzf ncbi-igblast-1.22.0-x64-linux.tar.gz
blastdbcmd -db ncbi-igblast-1.22.0/internal_data/human/human_TR_V -entry all > igblast.fa

mkdir refseq
cd refseq
for i in {1..14}; do
    wget -4 https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.${i}.protein.faa.gz
    wget -4 https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.${i}.rna.fna.gz
done
cat human.*.protein.faa.gz | gunzip > human.protein.faa
cat human.*.rna.fna.gz | gunzip > human.protein.fna
rm -f ./*.gz
cd ..

gunzip ./*.gz

wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz
tar xvzf refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz
python imgt_parse.py
python merge_data.py

for fn in out4msa.aa.fa.d/*.fa out4msa.nt.fa.d/*.fa; do
    mafft --thread -1 --auto "${fn}" > "${fn}".mafft
done

