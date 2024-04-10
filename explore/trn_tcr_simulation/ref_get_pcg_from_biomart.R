library("tidyverse")
library("arrow")
library("biomaRt")

ensembl97 <- useMart(
    host = 'https://jul2019.archive.ensembl.org',
    biomart = 'ENSEMBL_MART_ENSEMBL',
    dataset = 'hsapiens_gene_ensembl'
)
attrs <- listAttributes(ensembl97) %>% as_tibble()

gene_trans <- getBM(
    attributes = c("ensembl_gene_id_version", "ensembl_transcript_id_version"),
    mart = ensembl97
)
mane_trans <- getBM(
    attributes = c("transcript_mane_select", "ensembl_transcript_id_version"),
    mart = ensembl97
) %>%
    dplyr::filter(transcript_mane_select != "")
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
readr::write_tsv(selected_gene_trans, "ref/ens.trans_gene_map.tsv",col_names=FALSE)
