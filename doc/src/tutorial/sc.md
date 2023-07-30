# Tutorial on Simulation of scRNA-Seq Data

This is a tutorial on how to generate scRNA-Seq data from reference genome.

```{warning}
This version of scRNA-Seq simulator **does not** support:

- Generation of raw FASTQ reads from manufacturers like 10xGenomics, etc.
- Generation of Differentially Expressed Genes (DEGs)/Different Isoform Usage (DIU), or data distributed in Negative Biomial distribution as-is observed in real scRNA-Seq data.
```

## Environment Specifications

Here would list version information of each component used in this tutorial for reproductive purposes.

| Software     | Version   |
|--------------|-----------|
| [GNU Bash](https://www.gnu.org/software/bash)     | 5.2.15(1) |
| [GNU Grep](https://www.gnu.org/software/grep)     | 3.8       |
| [GNU Wget](https://www.gnu.org/software/wget)     | 1.21.3    |
| [yasim](https://pypi.org/project/yasim/)        | 3.1.5     |
| yasim\_sctcr | 0.1.0     |
| [art\_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) | 2.5.8     |
| [samtools](https://www.htslib.org)     | 1.17      |

Make sure you had correctly set up the environment as-is specified in {doc}`Readme </_root/Readme>`.

## Step 0. Preparation

Following code retrieves reference genome sequence (FASTA format) and annotations (GTF format) of human from UCSC used in this example. Only chromosome 14 and 7 are used to cut down simulation time.

```shell
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
# Filter chromosome 7 and 14 only
gzip -dcf hg38.ncbiRefSeq.gtf.gz \
    | grep -e '^chr7\s' -e '^chr14\s' \
    > hg38.ncbiRefSeq_chr7_14.gtf

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gzip -dcf hg38.fa.gz > hg38.fa
samtools faidx hg38.fa
```

## Step 1. Generate Cell Barcodes

Now we would try to generate barcodes for 10 cells. Barcode are unique identifiers for each cell.

```shell
python -m yasim_sc generate_barcode -n 10 -o barcode.txt
```

**Generates:**

- `barcode.txt`, which is a flat file of barcode nucleotide sequences separated by `\n`.

## Step 2. Sample Protein-Coding Genes

This step removes genes except protein-coding genes (PCGs) from reference genome and downsample number of PCGs to reduce simulation time. It is purely optional.

This step would require a `ncbi_dataset.tsv` available from [here](https://www.ncbi.nlm.nih.gov/labs/data-hub/gene/table/taxon/9606/). You can get version used in the article [here](../../../data/ncbi_dataset.tsv.xz).

The following example samples 20 protein coding genes from chromosome 7 and 14 for downstream simulation:

```shell
python -m yasim_sc sample_pcg \
    -i ncbi_dataset.tsv.xz \
    -g hg38.ncbiRefSeq_chr7_14.gtf \
    -o hg38.ncbiRefSeq_subsampled.gtf \
    --num_genes_to_sample 20
```

**Generates:**

- `hg38.ncbiRefSeq_subsampled.gtf`, the GTF containing selected genes.

## Step 3. Generate Gene-Level and Isoform-Level Depth

This step is similiar as-is in bulk version of YASIM. See documentation there for more details.

The following example generates gene- and isoform-level depth with mean sample-level sequencing depth 5.

```shell
python -m yasim generate_gene_depth \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -d 5 \
    -o notcr_gene_depth.tsv
python -m yasim generate_isoform_depth \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -d notcr_gene_depth.tsv \
    -o notcr_isoform_depth.tsv
```

**Generates:**

- `notcr_gene_depth.tsv`, which is a YASIM Gene-level-depth compatible depth file.
- `notcr_isoform_depth.tsv`, which is a YASIM Isoform-level-depth compatible depth file.

## Step 4. Generate Isoform-Level Depth Information for Each Cell

This step introduces cell-level differences by adding uniformally-distributed machine-noise.

```shell
python -m yasim_sc generate_barcoded_isoform_replicates \
    -d notcr_isoform_depth.tsv \
    -b tcr_barcode.txt \
    -o notcr_isoform_depth_sc.d
```

**Generates:**

- `notcr_isoform_depth_sc.d`, a directory where isoform-level depth of each cell are placed.

## Step 5. Transcribe into cDNA FASTA

Same as YASIM, so not introduced here:

```shell
python -m labw_utils.bioutils transcribe \
    -f hg38.fa \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -o notcr_trans.fa
```

## Step 6. Invocation of LLRGs

The following example would perform scRNA-Seq using ART.

```shell
python -m yasim_sc art \
    -F notcr_trans.fa.d \
    -o sim_notcr_50.d \
    -d notcr_isoform_depth_sc.d \
    -e art_illumina \
    -j 20
```

Generates:

- `sim_notcr_50.d`: The target directory where reads generated from all cells are separated into single file.
