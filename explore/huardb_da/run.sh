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
samtools faidx cdr3_aa.fa
seqtk sample all_contigs.fq 20000 >all_contigs_s20k.fq

mkdir -p ref
cd ref
axel https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
axel https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
axel https://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
singularity build deepvariant_cpu_1.5.0.sif docker://google/deepvariant:1.5.0

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

for fn in ens.cdna.fa.d/*.fa; do
    base_name="$(basename "${fn}")"
    echo "${base_name}"
    if [ ! -f ens.cdna.fa.d/"${base_name}.bwt" ]; then
        bwa index ens.cdna.fa.d/"$(basename "${fn}")"
    fi
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
    freebayes-parallel \
        <(fasta_generate_regions.py ens.cdna.fa.fai 1000) 40 \
        -f ens.cdna.fa \
        -p 1 \
        out4msa_samples.nt."${aligner}".rg.bam >out4msa_samples.nt."${aligner}".freebayes.vcf
    samtools mpileup out4msa_samples.nt."${aligner}".rg.bam -f ens.cdna.fa |
        varscan mpileup2cns \
            --output-vcf \
            --vcf-sample-list <(echo huardb) \
            --variants >out4msa_samples.nt."${aligner}".varscan.vcf
    # mamba create -n yasim-clair3 clair3 -c conda-forge -c bioconda
    run_clair3.sh \
        --bam_fn=out4msa_samples.nt."${aligner}".rg.bam \
        --ref_fn=ens.cdna.fa \
        --threads=40 \
        --include_all_ctgs \
        --platform=ilmn \
        --sample_name=huardb \
        --model_path="${CONDA_PREFIX}"/bin/models/ilmn/ \
        --output=out4msa_samples.nt."${aligner}".clair3.vcf.d
    zcat out4msa_samples.nt."${aligner}".clair3.vcf.d/merge_output.vcf.gz >out4msa_samples.nt."${aligner}".clair3.vcf
    # mamba create -n yasim-deepvariant deepvariant -c conda-forge -c bioconda etils python=3.7 # Failed

    # ERR: Some or all of the contigs in the reference genome (ens.cdna) are not present in the read files.
    # octopus --reference ens.cdna.fa --reads out4msa_samples.nt."${aligner}".rg.bam > out4msa_samples.nt."${aligner}".octopus.vcf
done
python collect_vcf.py >variants.tsv

for fn in out4msa_samples.nt.bwa.*.vcf; do
    bgzip "${fn}"
    tabix "${fn}".gz
done

mkdir -p out4msa_samples.nt.bwa.merged.d
bcftools isec \
    -p out4msa_samples.nt.bwa.merged.d \
    -n +2 \
    out4msa_samples.nt.bwa.*.vcf.gz >/dev/null

for fn in out4msa_samples.nt.bwa.merged.d/*.vcf; do
    bgzip "${fn}"
    tabix "${fn}".gz
done

bcftools concat -a --remove-duplicates out4msa_samples.nt.bwa.merged.d/*.vcf.gz >out4msa_samples.nt.bwa.merged.vcf

mkdir -p gene4msa
for gene_name in TRAV TRAJ TRBV TRBJ TRBC; do
    for type in cdna pep; do
        seqkit grep \
            --by-name \
            --use-regexp \
            -p "${gene_name}.*" \
            ens."${type}".fa \
            >gene4msa/"${gene_name}.${type}.fa"
        # t_coffee -type=dna -output fasta_aln gene4msa/"${gene_name}".fa >gene4msa/"${gene_name}.${type}".fa.t_coffee
        mafft --thread -1 --auto gene4msa/"${gene_name}.${type}".fa >gene4msa/"${gene_name}.${type}".mafft
        clustalo --auto -i gene4msa/"${gene_name}.${type}".fa >gene4msa/"${gene_name}.${type}".clustalo
        famsa gene4msa/"${gene_name}.${type}".fa gene4msa/"${gene_name}.${type}".famsa
        muscle -align gene4msa/"${gene_name}.${type}".fa -output gene4msa/"${gene_name}.${type}".muscle
        probcons gene4msa/"${gene_name}.${type}".fa >gene4msa/"${gene_name}.${type}".probcons
    done
done
mkdir -p mut_aa.d
