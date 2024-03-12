#!/usr/bin/env bash
set -ue
# On datacenter@10.106.125.76, run:
find /Volumes/rsch/HuarcBackup -name all_contig_annotations.json >huardb_all_contig_annotations.lofn

scp -P 10022 datacenter@10.106.125.76:/Volumes/rsch/HuarcBackup/huardb_all_contig_annotations.lofn .
mkdir -p raw_data

while read -r line; do
    tfn="raw_data/$(echo "${line}" | tr / _)"
    if [ ! -f "${tfn}" ]; then
        scp -P 10022 datacenter@10.106.125.76:"${line}" "raw_data/$(echo "${line}" | tr / _)"
    fi
done <huardb_all_contig_annotations.lofn

mkdir -p pq
python read_json.py
python merge_data.py

#for fn in out4msa_samples.nt.fa.d/*.fa; do
#    echo "${fn}"
#    if [ ! -f "${fn}".famsa ]; then
#        famsa -t 0 "${fn}" "${fn}".famsa
#    fi
#    hmmbuild --cpu 40 --fast --dna "${fn}".hmm "${fn}".famsa
#done
#for fn in out4msa_samples.nt.fa.d/*.hmm; do
#    GEN_URL="$(curl -H 'Accept:application/json' -F file="@${fn}" -F processing=hmm http://localhost:8000/ | jq .uuid --raw-output)"
#    echo "${GEN_URL}"
#    curl -H 'Accept:image/png' "http://localhost:8000/logo/${GEN_URL,,}/download/?colors=default&format=png" >"${fn}".png
#done
#
#Rscript plot.R
#
#bedtools sort -i all_contigs.cellranger.bed >all_contigs.cellranger.srt.bed
#
#seqtk sample all_contigs.fa 1000 >all_contigs.h1k.fa
#samtools faidx all_contigs.h1k.fa
#pblat -threads=36 -t=dna -q=dna -out=psl ../test/all_nt.fa all_contigs.h1k.fa all_contigs.h1k.nt.blat.fa.psl
#pslToBed all_contigs.h1k.nt.blat.fa.psl /dev/stdout | awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$7,$8,$1,".",$6}' >all_contigs.h1k.nt.blat.fa.bed
#mmseqs easy-search \
#    --format-output "query,target,evalue,cigar,alnlen,nident,qlen,tlen,qstart,qend,tstart,tend" \
#    -s 7.5 \
#    --search-type 4 \
#    --format-mode 0 \
#    all_contigs.h1k.fa ../test/all_nt.fa all_contigs.h1k.fa.mmseqs2.tsv /tmp
#cut all_contigs.h1k.fa.mmseqs2.tsv -f 1,2,9,10 |
#    awk 'BEGIN{FS="\t";OFS="\t"}{if($3>$4){print $1,$4,$3,$2,".","-";}else{print $1,$3,$4,$2,".","+";}}' \
#        >all_contigs.h1k.fa.mmseqs2.bed
#
#awk 'BEGIN{FS="\t";OFS="\t"}{print $1,"1",$2}' all_contigs.h1k.fa.fai >all_contigs.h1k.fa.fai.bed
#bedtools intersect \
#    -a all_contigs.cellranger.srt.bed \
#    -b all_contigs.h1k.fa.fai.bed -wa >all_contigs.h1k.cellranger.bed
#
#mmseqs easy-search \
#    --format-mode 1 \
#    -s 7.5 \
#    --search-type 3 \
#    out4msa_samples.nt.fa.d/TRAC.fa \
#    ../test/out4msa.nt.fa.d/TRAC.fa \
#    test.sam \
#    /tmp
#samtools sort test.sam -o test.bam -@40
#samtools index test.bam

mkdir -p ref
cd ref
axel https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
axel https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
axel https://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
gunzip ./*.gz

cd ..
for name in cdna pep; do
    seqkit grep \
        --by-name \
        --use-regexp \
        -p 'chromosome:GRCh38:(7|14):' \
        ref/Homo_sapiens.GRCh38."${name}".all.fa |
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
gatk CreateSequenceDictionary -R ens.cdna.fa
python -m labw_utils.bioutils split_fasta ens.cdna.fa

#gffread \
#    -F \
#    --bed \
#    --stream ref/Homo_sapiens.GRCh38.111.gtf |
#    grep -E 'gene_biotype=TR_[VDJC]_gene' |
#    grep -E 'gene_name=TR[AB]' |
#    grep -E '^([0-9]+|X|Y|Mt)\s' |
#    cut -f 1,6,7,8,13 |
#    sed -E 's/CDS.*gene_name=([^;]+);.*/\1/' \
#        >ens.tcrs.bed
#printf >ens.tcrs_chr.bed
#
#while IFS=$'\t' read -a line; do
#    chr="${line[0]}"
#    chr_start="${line[2]}"
#    chr_end="${line[3]}"
#    gene_name="${line[4]}"
#    if [ ! -f "${gene_name}".dbsnp.ljson ]; then
#        echo "${gene_name} -- ${chr}:${chr_start}-${chr_end}"
#        esearch -db snp -query "${chr}"'[Chromosome] AND ('"${chr_start}"':'"${chr_end}"'[Base Position])' | efetch -format json >"${gene_name}".dbsnp.ljson
#    fi
#done <ens.tcrs.bed
# Incorrect: introns included as well.

for fn in ens.cdna.fa.d/*.fa; do
    base_name="$(basename "${fn}")"
    echo "${base_name}"
    if [ ! -f ens.cdna.fa.d/"${base_name}.bwt" ]; then
        bwa index ens.cdna.fa.d/"$(basename "${fn}")"
    fi
    #    if [ -f out4msa_samples.nt.fa.d/"${base_name}" ] && [ ! -f out4msa_samples.nt.fa.d/"${base_name}".mmseqs2.bam ]; then
    #        mmseqs easy-search \
    #            --format-mode 1 \
    #            --search-type 3 \
    #            out4msa_samples.nt.fa.d/"${base_name}" \
    #            ens.cdna.fa.d/"${base_name}" \
    #            out4msa_samples.nt.fa.d/"${base_name}".mmseqs2.sam \
    #            /tmp
    #        samtools sort \
    #            -@40 out4msa_samples.nt.fa.d/"${base_name}".mmseqs2.sam \
    #            -o out4msa_samples.nt.fa.d/"${base_name}".mmseqs2.bam
    #        samtools index out4msa_samples.nt.fa.d/"${base_name}".mmseqs2.bam
    #    fi
    base_name_fq="${base_name/fa/fq}"
    echo "${base_name_fq}"
    if [ -f out4msa_samples.nt.fq.d/"${base_name_fq}" ] &&
        [ "$(wc -l out4msa_samples.nt.fq.d/"${base_name_fq}" | cut -f 1 -d ' ')" -gt 400 ] &&
        [ ! -f out4msa_samples.nt.fq.d/"${base_name_fq}".bwa.bam ]; then
        bwa mem -x intractg -t 40 ens.cdna.fa.d/"${base_name}" out4msa_samples.nt.fq.d/"${base_name_fq}" >out4msa_samples.nt.fq.d/"${base_name_fq}".bwa.sam
        samtools sort \
            -@40 out4msa_samples.nt.fq.d/"${base_name_fq}".bwa.sam \
            -o out4msa_samples.nt.fq.d/"${base_name_fq}".bwa.bam
        samtools index out4msa_samples.nt.fq.d/"${base_name_fq}".bwa.bam
    fi
done

# samtools merge -f -@40 out4msa_samples.nt.mmseqs2.bam out4msa_samples.nt.fa.d/*.mmseqs2.bam
samtools merge -f -@40 out4msa_samples.nt.bwa.bam out4msa_samples.nt.fq.d/*.bwa.bam

for aligner in bwa; do
    samtools index out4msa_samples.nt."${aligner}".bam
    samtools addreplacerg \
        -r 'PL:10x_genomics' \
        -r 'ID:0' \
        -r "LB:huardb.${aligner}" \
        -r 'SM:huardb' \
        -r 'PU:na' \
        -o out4msa_samples.nt."${aligner}".rg.bam \
        out4msa_samples.nt."${aligner}".bam
    samtools index out4msa_samples.nt."${aligner}".rg.bam
    gatk HaplotypeCaller \
        -R ens.cdna.fa \
        -ploidy 1 \
        -I out4msa_samples.nt."${aligner}".rg.bam \
        -O out4msa_samples.nt."${aligner}".gatk.vcf
    bcftools mpileup -Ou \
        out4msa_samples.nt."${aligner}".rg.bam \
        -f ens.cdna.fa |
        bcftools call \
            --ploidy 1 \
            --multiallelic-caller \
            --variants-only \
            --threads 40 \
            -o out4msa_samples.nt."${aligner}".bcftools.vcf
    freebayes-parallel \
        <(fasta_generate_regions.py ens.cdna.fa.fai 1000) 40 \
        -f ens.cdna.fa \
        -p 1 \
        out4msa_samples.nt."${aligner}".rg.bam >out4msa_samples.nt."${aligner}".freebayes.vcf
done
