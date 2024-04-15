library(tidyverse)

df <- readr::read_tsv(
    "pbmc_1k_v3_possorted_genome_bam.salmon.sorted.depth.tsv.xz",
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
        mean_depth>200,
        len>3000
    ) %>%
    dplyr::sample_n(20)

df_sel <- df %>%
    dplyr::filter(CHROM%in%df_mean_ex$CHROM) %>%
    dplyr::group_by(CHROM)


ggplot(df_sel) +
    geom_line(aes(
        x=POS,
        y=DEPTH
    )) +
    facet_wrap(.~CHROM, scales = "free") +
    theme_bw()

