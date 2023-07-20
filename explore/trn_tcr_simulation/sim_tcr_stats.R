library(tidyverse)
library(pheatmap)
library(ggalluvial)
df <- readr::read_tsv("sim_tcr.stats.tsv", quote = "'")

dfm <- df %>%
    dplyr::mutate(
        ACDR3_AA_LEN = sapply(as.vector(.$ACDR3_AA), stringr::str_length),
        BCDR3_AA_LEN = sapply(as.vector(.$BCDR3_AA), stringr::str_length),
        ACDR3_AA_TERM = as.vector(sapply(
            as.vector(.$ACDR3_AA),
            function(x) {
                stringr::str_sub(x, start = stringr::str_length(x))
            }
        ))
    )
ggplot(dfm) + geom_histogram(aes(x = BCDR3_AA_LEN), binwidth = 1)
ggplot(dfm) + geom_histogram(aes(x = ACDR3_AA_LEN), binwidth = 1)

freq_df <- dfm %>%
    dplyr::group_by(
        ACDR3_AA_TERM, TRAJ
    ) %>%
    dplyr::summarise(n = n()) %>%
    tidyr::pivot_wider(
        id_cols = TRAJ,
        names_from = ACDR3_AA_TERM,
        values_fill = 0,
        values_from = n
    ) %>%
    as.data.frame()
rownames(freq_df) <- freq_df$TRAJ
freq_df$TRAJ <- NULL
pheatmap(freq_df)
