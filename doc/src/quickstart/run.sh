#!/usr/bin/env bash
set -ue

if [ ! -f ce11.ncbiRefSeq.chr1.gtf ]; then
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz &>>/dev/null
    gunzip -f ce11.ncbiRefSeq.gtf.gz
    grep -i '^chrI\s' <ce11.ncbiRefSeq.gtf >ce11.ncbiRefSeq.chr1.gtf
fi
if [ ! -f ce11.chr1.fa ]; then
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz &>>/dev/null
    gunzip -f ce11.fa.gz
    head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1) - 1))" >ce11.chr1.fa
fi
function grep_pbar() {
    { grep -v '1%' || true; } |
        { grep -v '26%' || true; } |
        { grep -v '51%' || true; } |
        { grep -v '76%' || true; }
}
{
    if [ ! -f ce11.ncbiRefSeq.chr1.as.gtf ]; then
        python -m yasim generate_as_events \
            -f ce11.fa \
            -g ce11.ncbiRefSeq.chr1.gtf \
            -o ce11.ncbiRefSeq.chr1.as.gtf \
            -c 5 2>&1
    fi
} | {
    grep -v 'inferred from feature transcript' || true
} | grep_pbar >>generate_as_events.log
{
    if [ ! -f ce11_gene_depth.tsv ]; then
        python -m yasim generate_gene_depth \
            -g ce11.ncbiRefSeq.chr1.as.gtf \
            -o ce11_gene_depth.tsv \
            -d 5 2>&1
    fi
} | grep_pbar &>>generate_gene_depth.log
{
    if [ ! -f ce11_isoform_depth.tsv ]; then
        python -m yasim generate_isoform_depth \
            -g ce11.ncbiRefSeq.chr1.as.gtf \
            -d ce11_gene_depth.tsv \
            -o ce11_isoform_depth.tsv \
            --alpha 4 2>&1
    fi
} | grep_pbar &>>generate_isoform_depth.log
{
    if [ ! -f ce11_trans_as.fa ]; then
        python -m labw_utils.bioutils transcribe \
            -f ce11.chr1.fa \
            -g ce11.ncbiRefSeq.chr1.as.gtf \
            -o ce11_trans_as.fa 2>&1
    fi
} | grep_pbar &>>transcribe.log
{
    if [ ! -f pbsim3_mode.fq ]; then
        python -m yasim pbsim3 \
            -F ce11_trans_as.fa.d \
            -j 40 \
            -e /home/yuzj/bin/pbsim3 \
            --ccs_pass 20 \
            -d ce11_isoform_depth.tsv \
            -m RSII \
            -M qshmm \
            --strategy trans \
            -o pbsim3_mode \
            2>&1
    fi
} | {
    grep -v WARNING || true
} | grep_pbar &>>pbsim3.log
