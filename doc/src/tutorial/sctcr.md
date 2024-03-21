# Tutorial on Simulation of scTCR-Seq Data

This is a tutorial on how to generate scTCR-Seq data from reference genome.

```{warning}
This version of scTCR-Seq simulator **does not** support clonal expansion or generation of raw FASTQ reads from manufacturers like 10xGenomics, etc.
```

Download of real statistical data:

- [cdr3_deletion_table.min.json.xz](../../../data/cdr3_deletion_table.min.json.xz): Human CDR3 deletion table.
- [cdr3_insertion_table.min.json.xz](../../../data/cdr3_insertion_table.min.json.xz): Human CDR3 insertion table.
- [tcr_cache.min.json.xz](../../../data/tcr_cache.min.json.xz): TCR Cache built against [UCSC hg38](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz) reference genome sequence and [UCSC NCBIRefSeq](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz) reference genome annotation.
- [tcr_genelist.min.json.xz](../../../data/tcr_genelist.min.json.xz): Human TCR gene list curated from IMGT.
- [IMGT_Protein_Display.min.json.xz](../../../data/IMGT_Protein_Display.min.json.xz): Human TCR reference Amino-Acid (AA) Sequence curated from IMGT.

## Environment Specifications

Here would list version information of each component used in this tutorial for reproductive purposes.

| Software                                                                                           | Version   |
|----------------------------------------------------------------------------------------------------|-----------|
| [GNU Bash](https://www.gnu.org/software/bash)                                                      | 5.2.15(1) |
| [GNU Grep](https://www.gnu.org/software/grep)                                                      | 3.8       |
| [GNU Wget](https://www.gnu.org/software/wget)                                                      | 1.21.3    |
| [yasim](https://pypi.org/project/yasim/)                                                           | 3.1.5     |
| [yasim\_sctcr](https://pypi.org/project/yasim-sctcr/)                                              | 0.1.0     |
| [art\_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) | 2.5.8     |
| [samtools](https://www.htslib.org)                                                                 | 1.17      |
| [seqkit](https://bioinf.shenwei.me/seqkit)                                                         | XXX       |

Make sure you had correctly set up the environment as-is specified in {doc}`Readme </_root/Readme>`.

## Step 0. Preparation

The following code generates FASTA for peptides and cDNA sequences for TRAV/TRBV genes used in this simulator from Ensembl. The simulator uses data built from CellRanger `vdj` using Ensembl 97 for reference, so we use this version in our tutorial. You may use later versions instead.

```shell
wget https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
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
```

The above command will generate `ens.cdna.fa` for cDNA sequences and `ens.pep.fa` for peptide sequences and their indices in FAI format if no error occurs.

## Step 1. Generate Cell Barcodes

Now we would try to generate barcodes for 10 cells. Barcodes are unique identifiers for each cell.

```shell
python -m yasim_sc generate_barcode -n 10 -o barcode.txt
```

**Generates:**

- `barcode.txt`, which is a flat file of barcode nucleotide sequences separated by `\n`.

## Step 2. Generate TCR Depth

The Following code generates a depth file with 20 depth for each cell.

```shell
python -m yasim_sctcr generate_tcr_depth \
    -b barcode.txt \
    -o tcr_depth.tsv \
    -d 20
```

**Generates:**

- `tcr_depth.tsv`, which is a YASIM Isoform-level-depth compatible depth file. Here, TRA and TRB contig of a cell is regarded as one isoform.

## Step 3. Generate TCR Cache

The TCR cache is a JSON containing TCR Nucleotide (NT) to Amino-Acid (AA) alignment information, which is required for simulation of TCR rearrangement. The generated file can be reused within same species.

```shell
python -m yasim_sctcr generate_tcr_cache \
    --tcr_aa_table_path IMGT_Protein_Display.min.json.xz \
    --tcr_genelist_path tcr_genelist.min.json.xz \
    -f hg38.fa \
    -g hg38.ncbiRefSeq_chr7_14.gtf \
    -o tcr_cache.json
```

This generates following files:

- `tcr_cache.json`, the generated TCR cache. User may not use this file.
- `tcr_cache.json.aa.fa` and `tcr_cache.json.nt.fa`, NT and AA sequence of generated TCR cache. For debug purpose only;User may not use this file.

## Step 4. Simulate Ground-Truth TCR Contigs

This step may report `finish with N failures` -- don't worry! None of your cells would lost.

```shell
python -m yasim_sctcr rearrange_tcr \
    --tcr_genelist_path tcr_genelist.min.json.xz \
    --cdr3_deletion_table_path cdr3_deletion_table.min.json.xz \
    --cdr3_insertion_table_path cdr3_insertion_table.min.json.xz \
    --tcr_cache_path tcr_cache.json \
    -b barcode.txt \
    -o sim_tcr
```

This generates following files:

- `sim_tcr.stats.tsv`, the ground-truth TCR contig statistics, with following columns:
  - `UUID`, the barcode, which is previously in UUID format.
  - `TRAV` \& `TRAJ` \& `TRBV` \& `TRBJ`, gene name of corresponding TCR segments.
  - `ACDR3_AA` \& `BCDR3_AA`, AA sequence of CDR3 on corresponding TCR segments.
  - `ACDR3_NT` \& `BCDR3_NT`, NT sequence of CDR3 on corresponding TCR segments.
  - `ALPHA_AA` \& `BETA_AA`, AA sequence of corresponding TCR chain.
  - `ALPHA_NT` \& `BETA_NT`, NT sequence of corresponding TCR chain.
- `sim_tcr.aa.fa` and `sim_tcr.nt.fa`, ground-truth contig in AA and NT with seqname `{barcode}_A` for TCR alpha chain and `{barcode}_B` for TCR beta chain.

## Step 5. Invocation of LLRGs

The following example would perform single-end scTCR-Seq using ART.

```shell
# Split the FASTA before performing calling YASIM RNA-Seq interface
python -m labw_utils.bioutils split_fasta sim_tcr.nt.fa
python -m yasim art \
    -F sim_tcr.nt.fa.d \
    -o sim_tcr_50 \
    --sequencer_name GA2 \
    --read_length 50 \
    -d tcr_depth.tsv \
    -e art_illumina \
    --preserve_intermediate_files \
    -j 20
```

Generates:

- `sim_tcr_50.d`: The target directory where reads generated from TRA and TRB sequences of all cells are separated into single file.
- `sim_tcr_50.fq`: A merged file which contains all simulated reads.
- `sim_tcr_50.fq.stats`: Real simulation statistics.
