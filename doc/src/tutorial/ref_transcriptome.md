# Reproducibility Guide on Building References

This guide allows you to recreate the human reference published at our [Zenodo](https://doi.org/10.5281/zenodo.12155540) repository.

```{note}
Modern GNU/Linux is assumed for this tutorial. If you're working under Microsoft Windows, you may use [WSL](https://learn.microsoft.com/en-us/windows/wsl/).

Basic knowledge on Shell programming, R (especially [Tidyverse](https://www.tidyverse.org/) series and [Seurat](https://satijalab.org/seurat/)) and Python are assumed for this tutorial.

You should also know the basis of [regular expressions](https://en.wikipedia.org/wiki/Regular_expression) and its [Python implementation](https://docs.python.org/3/library/re.html).
```

Before the analysis, download (using GNU WGet for example) and extract the following files:

```shell
wget https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
```

## Making the Reference Transcriptome

We already know that YASIM LLRGs works on transcript-level, while scRNA-Seq data are usually provided on a gene-level basis. So, it is important to convert gene IDs to its **representative** transcripts. Here we will use protein-coding genes (PCGs) in Ensembl release 97 that have HGNC symbol and use MANE-selected transcripts as representative transcripts. If a gene has multiple MANE-selected transcripts, we will randomly decide one.

Firstly, we need a mapping between HGNC gene symbols and transcript IDs. The following R script could perform this:

```r
library("tidyverse")
library("biomaRt")

# You may use other Ensembl release version if you like.
ensembl97 <- useMart(
    host = 'https://jul2019.archive.ensembl.org',
    biomart = 'ENSEMBL_MART_ENSEMBL',
    dataset = 'hsapiens_gene_ensembl'
)

# Gene ID-Transcript ID mapping.
gene_trans <- getBM(
    attributes = c("ensembl_gene_id_version", "ensembl_transcript_id_version"),
    mart = ensembl97
)
# Ensembl transcript ID-MANE selected RefSeq ID mapping.
mane_trans <- getBM(
    attributes = c("transcript_mane_select", "ensembl_transcript_id_version"),
    mart = ensembl97
) %>%
    dplyr::filter(transcript_mane_select != "")
# Ensembl gene ID-HGNC symbol mapping.
hgnc_gene <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id_version", "gene_biotype"),
    mart = ensembl97
) %>%
    dplyr::filter(hgnc_symbol != "", gene_biotype == "protein_coding")

selected_gene_trans <- gene_trans %>%
    dplyr::inner_join(mane_trans, by="ensembl_transcript_id_version") %>%
    dplyr::select(!transcript_mane_select) %>%
    dplyr::inner_join(hgnc_gene, by="ensembl_gene_id_version") %>%
    dplyr::select(!c(ensembl_gene_id_version, gene_biotype))
readr::write_tsv(selected_gene_trans, "ens.trans_gene_map.tsv", col_names=FALSE)
```

**Generates:**

- `ens.trans_gene_map.tsv`, which is a Tab-Separated Value (TSV) whose first column is Ensembl transcript ID and second column is HGNC symbol of corresponding gene.

Now we would download Ensembl transcript cDNAs and subset them to exclude unselected transcripts.

```shell
# Use sektk to select selected transcripts.
seqtk subseq \
    Homo_sapiens.GRCh38.cdna.all.fa \
    <(cut -f 1 ens.trans_gene_map.tsv) \
    >ens.sel_genes.cdna.fa
```

**Generates:**

- `ens.sel_genes.cdna.fa`, a FASTA of all selected transcripts that is considerably smaller than the original one.

## Constructing the TCR Cache

The TCR cache is a **reusable** JSON containing TCR nucleotide to amino acid alignment information, which is required for simulation of TCR rearrangement. The generated file can be reused for all studies using the same reference TCR sequence.

Before generating the TCR cache, we would download TCR sequences in cDNA and amino acid from Ensembl. The following code generates FASTA for peptides and cDNA sequences for TRAV/TRBV genes used in this simulator from Ensembl. The simulator uses data built from huARdb v1 which uses Ensembl 97 for reference, so we use this version in our tutorial.

```shell
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
            /dev/stdin | 
        sed -E 's;^>.+ gene_symbol:(\S+) .+$;>\1;' \
        >ens."${name}".fa
done
```

The above command will select genes that:

- Are on chromosome 7 or 14, which excludes pesudogene `TRBV20OR9-2` (see it on [ENSG00000205274](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000205274;r=9:33617845-33618508;t=ENST00000379435;redirect=no) or [HGNC:12197](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:12197)).
- Have `TR_V_gene`/`TR_D_gene`/`TR_J_gene`/`TR_C_gene` as their [biotype](https://www.ensembl.org/info/genome/genebuild/biotypes.html).
- Gene symbol matches `TR[AB]` regular expression.

Since there will be only 1 transcript for all TCR gene segments in Ensembl release 97, we replaced FASTA header for each record with the corresponding gene symbol. **This may not be true in other reference genomes or further releases of Esnembl.**

It will finally generate `ens.cdna.fa` for cDNA sequences and `ens.pep.fa` for peptide sequences if no error occurs. We may generate a TCR cache using:

```shell
python -m yasim_sctcr generate_tcr_cache \
    --tcr_cdna_fa_path ens.cdna.fa \
    --tcr_pep_fa_path ens.pep.fa \
    -o tcr_cache.json
```

**Generating:**

- `tcr_cache.json`, the generated TCR cache. User may not use this file.
