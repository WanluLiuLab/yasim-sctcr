library(tidyverse)

df <- readr::read_tsv(
    "SC5pv2_GEX_Human_Lung_Carcinoma_DTC_possorted_genome_bam.salmon.sorted.depth.tsv.xz",
    col_names = c("CHROM", "POS", "DEPTH"),
    skip=1
)

df_mean_ex <- df %>%
    dplyr::group_by(CHROM) %>%
    dplyr::summarise(
        mean_depth=mean(DEPTH),
        len=max(POS)
    ) %>%
    dplyr::filter(
        mean_depth>100,
        len>3000
    )

df_sel <- df %>%
    dplyr::filter(CHROM%in%df_mean_ex$CHROM) %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(POS=POS/max(POS), MAXDEPTH=max(DEPTH)) %>%
    dplyr::filter(DEPTH==MAXDEPTH)

ggplot(df_sel) +
    geom_histogram(aes(
        x=POS
    )) +
    theme_bw()

