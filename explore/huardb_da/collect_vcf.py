import sys

import pysam

if __name__ == "__main__":
    print("caller", "CHR", "START", "STOP", "REF", "ALT", "REF_COUNT", "ALT_COUNT", sep="\t")

    for caller in ["clair3", "freebayes", "gatk", "varscan"]:
        file_name = f"out4msa_samples.nt.bwa.{caller}.vcf"
        with pysam.VariantFile(file_name) as varr:
            for record in varr:
                if record.filter:
                    if len(record.filter) == 1 and record.filter[0].name == "PASS":
                        pass
                    else:
                        continue
                if record.samples[0].get("DP") < 100:
                    continue
                if len(record.alleles) != 2:
                    continue
                if len(record.samples) != 1:
                    continue
                if caller == "varscan":
                    ad = (record.samples[0].get("RD"), record.samples[0].get("AD"))
                else:
                    ad = record.samples[0].get("AD")
                print(caller, record.contig, record.start, record.stop, *record.alleles, *ad, sep="\t")
