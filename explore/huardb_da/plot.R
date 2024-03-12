library(tidyverse)
library(patchwork)

df <- arrow::read_parquet("merged.parquet") %>%
    tibble::as_tibble() %>%
    dplyr::select(!c("nt", "quals"))

df_v_grouped <- df %>% dplyr::select(sample, productive) %>%
    dplyr::group_by(sample, productive) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = productive, values_from = n) %>%
    dplyr::filter((`FALSE`+`TRUE`) > 100) %>%
    dplyr::transmute(sample, prod_rate=1.0 * `TRUE`/(`FALSE`+`TRUE`))

ggplot(df_v_grouped) +
    geom_violin(aes(x='1', y=prod_rate)) +
    theme_bw()

ggplot(df) +
    geom_bar(aes(x=productive)) +
    facet_wrap(.~v, scales = "free") +
    theme_bw()

df_mut <- df %>%
    dplyr::mutate(vjlen = as.integer(j_start - v_end)) %>%
    tidyr::replace_na(list(d="NULL"))
ggplot(df_mut) +
    geom_violin(aes(x=vjlen, y=productive), kernel = "gaussian", adjust=4) +
    theme_bw()


# assemb_full_df <- data.frame()
# for (fn in c("all_long_d.aa.fa.tsv", "all_long_d.nt.fa.tsv")) {
#     message(fn)
#     df <- readr::read_tsv(
#         fn, quote = "'", show_col_types = FALSE,
#         col_names = c("query", "target", "evalue", "cigar", "alnlen", "nident", "qlen", "tlen", "qstart", "qend", "tstart", "tend")
#     )
#     if (length(grep("aa", fn)) != 0) {
#         df <- df %>%
#             dplyr::mutate(tlen=tlen*3, alnlen=alnlen*3)
#     }
#     df <- df %>%
#         dplyr::mutate(
#             aln_q = alnlen / qlen,
#             aln_t = alnlen / tlen
#         ) %>%
#         dplyr::filter(evalue<1E-5)
#
#     assemb_full_df <- rbind(
#         assemb_full_df,
#         df %>% dplyr::mutate(fn = fn)
#     )
# }
# g1 <- ggplot(assemb_full_df) +
#     geom_hex(
#         aes(x = aln_q, y = aln_t)
#     ) +
#     xlab("Portion of Query") +
#     facet_grid(fn ~ .) +
#     scale_fill_continuous(trans = "log10") +
#     theme_bw() +
#     scale_y_continuous("Portion of Target", breaks = scales::breaks_extended(20), limits = c(0, 1)) +
#     scale_x_continuous("Portion of Query", breaks = scales::breaks_extended(20), limits = c(0, 1))
# g2 <- ggplot(assemb_full_df) +
#     geom_violin(aes(y = aln_q, x = fn)) +
#     theme_bw() +
#     scale_y_continuous("Portion of Query", breaks = scales::breaks_extended(20), limits = c(0, 1))
# g3 <- ggplot(assemb_full_df) +
#     geom_violin(aes(y = aln_t, x = fn)) +
#     theme_bw() +
#     scale_y_continuous("Portion of Target", breaks = scales::breaks_extended(20), limits = c(0, 1))
# g4 <- ggplot(assemb_full_df) +
#     geom_hex(aes(y = alnlen, x = qlen)) +
#     facet_grid(. ~ fn) +
#     scale_fill_continuous(trans = "log10") +
#     theme_bw() +
#     ggtitle("Portion of Query")
# g <- (g1 + (g2 / g3)) / g4
# ggsave(paste0("aligned_samples.pdf"), g, width = 10, height = 15)
#

