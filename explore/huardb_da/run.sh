#!/usr/bin/env bash
# On datacenter@10.106.125.76, run:
find /Volumes/rsch/HuarcBackup -name all_contig_annotations.json > huardb_all_contig_annotations.lofn

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

for fn in out4msa_samples.nt.fa.d/*.fa; do
    echo "${fn}"
    if [ ! -f "${fn}".famsa ]; then
        famsa -t 0 "${fn}" "${fn}".famsa
    fi
    hmmbuild --cpu 40 --fast --dna "${fn}".hmm "${fn}".famsa
done

