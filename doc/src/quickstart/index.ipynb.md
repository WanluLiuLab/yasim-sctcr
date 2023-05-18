---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# YASIM Quickstart

This is a quickstart documentation for YASIM. In this documentation, you would generate Third-generation Sequencing (TGS) RNA-Seq reads from _C. Elegans_ (worm) reference genome with Alternative Splicing (AS) events.

This tutorial assumes basic understandings on RNA-Seq and Shell scripting. Those `# SKIP` on each cell is used for acceleration, so please do not copy-and-ans-paste them into your terminal.

+++

**How to read this documentation**:

Code block with leading `%%bash` are Shell code blocks. For example:

```{code-cell}
%%bash
ls -lFh | grep ipynb
```

## YASIM Data Flow Diagram

The following diagram lists full YASIM workflow. You are not limited to this and can start at any position.

```{figure} ../fig/yasim_toplevel.svg
:width: 100%
:align: left
:alt: Common YASIM Workflow
```

+++

## Environment Specifications

Here would list version information of each component used in this tutorial for reproductive purposes.

```{note}
we would assume that PBSIM3 is installed at `/home/yuzj/bin/pbsim3`. Change that to your own path on execution.
```

TODO: To be changed

```{code-cell}
%%bash
bash --version | head -n 1
grep --version | head -n 1
axel --version | head -n 1

python -m yasim --version 2> /dev/null
```

## Step 0. Retrive Reference Genome Sequence and Annotation

Inside the example, chromosome 1 of _C. Elegans_ reference genome sequence (in FASTA format) and annotation (in GTF format) from UCSC is used.

Get first chromosome of CE11 reference genome sequence and annotation from UCSC.

```{code-cell}
:tags: [skip-execution]

%%bash
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
gunzip ce11.ncbiRefSeq.gtf.gz
grep -i '^chrI\s' < ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr1.gtf

axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ce11.fa.gz
head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" >  ce11.chr1.fa
```

## Step 1. Generation of AS Events: `generate_as_events`

This step would generate alternative splicing events. It would take reference genome annotation as input and generate GTF with AS events as output.

The following example generates AS events from chromosome 1 of CE11 reference genome annotation with complexity index 5 (Some warnings generated during GTF parsing and process bar was filtered out):

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim generate_as_events \
    -f ce11.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11.ncbiRefSeq.chr1.as.gtf \
    -c 5
```

Generates:

- `ce11.ncbiRefSeq.chr1.as.gtf`, the generated GTF with AS events.
- `ce11.ncbiRefSeq.chr1.gtf.0.4.gvpkl.xz`, if not exist. This is a cache file for the `labw_utils` GTF parser.

The generated GTF (V2API) or part of generated GTF which expresses (V3API) should be seen as ground truth for benchmarking AS detectors. If you wish to benchmark quantifiers only using V3API, this step can be safely omitted. The reference genome should have a complexity index between 1 and 2 with real organism in about 2.

+++

## Step 2. Generate Sequencing Depth of Gene: `generate_gene_depth`

This step would generate base gene expression level in coverage for each gene over some GTF.

```{warning}
The generated coverage is **NOT** number of reads generated! It cannot be used as ground truth to assess quantification software! The number of reads ground truth will be provided by LLRG UIs introduced below.
```

Example of coverage generation on genome annotation with AS events just generated with mean depth 5:

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim generate_gene_depth \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -o ce11_gene_depth.tsv \
    -d 5
```

Generates:

- `ce11_gene_depth.tsv`, a TSV file with following columns:
  - `GENE_ID`, the `gene_id` field in GTF.
  - `DEPTH`, gene expression level in coverage.
- `ce11.ncbiRefSeq.chr1.as.gtf.0.4.gvpkl.xz`, if not present.

+++

## Step 3. Generate Sequencing Depth of Isoform: `generate_isoform_depth`

Here assigns expression level in coverage to isoforms. The mean expression level of each isoform in some gene should be equal to the base expression level assigned to that gene in previous step.

Following is a typical example. The parameter `alpha` (default to 4) was used to adjust evenness between isoform expression levels inside a gene.

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim generate_isoform_depth \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -d ce11_gene_depth.tsv \
    -o ce11_isoform_depth.tsv
```

Generates:

- `ce11_isoform_depth.tsv`, a TSV file with following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `DEPTH`, isoform expression level in coverage.

+++

## Step 4. Transcribe GTF to FASTA: `transcribe`

This step is general-purpose. It can be used to transcribe (**stranded**) any GTF that contains some isoforms into a FASTA file of all cDNA and a directory with all cDNAs in separate FASTA files, and skip those isoforms whose region is not defined in genomic sequence FASTA file. It would not add post-transcriptional modifications.

For those who is familiar with [BedTools](https://bedtools.readthedocs.io), It should generate similar output with:

```shell
bedtools getfasta -nameOnly -s -fi [FASTA] -bed [GTF] > [OUT]
```

Example:

```{code-cell}
:tags: [skip-execution]

%%bash
python -m labw_utils.bioutils transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -o ce11_trans_as.fa
```

Generates:

- `ce11_transcripts.fa`, the generated cDNA sequence FASTA.
- `ce11_transcripts.fa.d`, the directory where every cDNA is stored as separate FASTA.
- `ce11_transcripts.fa.stats`, a TSV file with following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `GENE_ID`, the `gene_id` field in GTF.
  - `SEQNAME`, chromosome or contig name.
  - `START`, the `start` field in GTF, 1-based inclusive.
  - `END`, the `end` field in GTF, 1-based inclusive.
  - `STRAND`, the `strand` field in GTF.
  - `ABSOLUTE_LENGTH`, is `START` - `END` + 1.
  - `TRANSCRIBED_LENGTH`, length of the cDNA without introns and UTRs.
  - `GC`, GC content of the cDNA in percentage.

+++

## Step 5. Invocation of LLRGs

The Low-Level Read Generators (LLRGs) are programs that simulates DNA-Seq on some reference genome sequence by clipping reads in appropriate length and introducing sequencing errors. YASIM invokes LLRG on stranded cDNA sequences to generate RNA-Seq data. Here we would demonstrate their usage with following examples:

+++

### Invocation of TGS LLRGs: Use `pbsim3` for Example

```{warning}
The official build of PBSIM, PBSIM2 and PBSIM3 shares a common executable anme (`pbsim`) but with different argument layout. For convenience, I renamed executable of PBSIM2 to `pbsim2` and PBSIM3 to `pbsim3`. If you do not use this in your computer, please use the `-e` option.
```

`pbsim3` is a general-purposed TGS DNA- and RNA-Seq simulator that supports multiple PacBio and Oxford Nanopore sequencers. It can generate Circular Consensus Sequence (CCS)/HiFi data.

Compared to NGS simulators, TGS simulators have `truncate_ratio_3p` and `truncate_ratio_5p`. These two parameters are used to set hard limits at two sides that allows simulation of incomplete reads due to reasons like 3' truncation.

Following is an example of simulation of CCS data using PacBio RS II model:

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim pbsim3 \
    -F ce11_trans_as.fa.d \
    -j 40 \
    -e /home/yuzj/bin/pbsim3 \
    --ccs_pass 20 \
    -d ce11_isoform_depth.tsv \
    -m RSII \
    -M qshmm \
    --strategy trans \
    -o pbsim3_mode
```

Generates:

- `pbsim3_mode.fq`, simulated Single-End FASTQ.
- `pbsim3_mode.d`, temporary directory that can be safely deleted.
- `pbsim3_mode.fq.stats`, statistics of simulated FASTQ. a TSV containing following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `INPUT_DEPTH`, isoform expression level in coverage provided by the upstream source.
  - `SIMULATED_N_OF_READS`, simulated number of reads. **This value can be used in assessing read quantifiers**.
  - `SIMULATED_N_OF_BASES`, simulated number of bases.
  - `TRANSCRIBED_LENGTH`, length of the cDNA without introns and UTRs.
  - `SIMULATED_DEPTH`, simulated depth.
