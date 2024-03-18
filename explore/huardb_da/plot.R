library(tidyverse)
library(patchwork)
library(pheatmap)
library(ComplexUpset)

df <- arrow::read_parquet("merged.parquet") %>%
    tibble::as_tibble() %>%
    dplyr::select(!c("nt", "quals"))

df_mut <- readr::read_tsv("/home/yuzj/Documents/tcrAnnotator/tmp/mut.tsv")

for (gene_name in unique(df_mut$gene)){
    message(gene_name)
    if (length(grep("CDR3", gene_name)) != 0) {
        next
    }
    p <- df_mut %>%
        dplyr::filter(gene==gene_name) %>%
    ggplot() +
    geom_bar(aes(fill=mutType, x=pepPos)) +
    theme_bw()
    ggsave(sprintf("mut_aa.d/%s.png", gene_name), p, limitsize = FALSE)
}

quit()


df_frame <- readr::read_tsv("/home/yuzj/Documents/tcrAnnotator/tmp/frame.csv")


df_pq <- arrow::read_parquet("merged_pq_cell.parquet") %>%
    tibble::as_tibble() %>%
    dplyr::select(!"sample_level_barcode")
df_clonal <- df_pq %>%
    dplyr::group_by(nt_blake2b) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n!=1) %>%
    dplyr::arrange(n) %>%
    dplyr::mutate(rank=n() - 1:n()) %>%
    tibble::as_tibble()
# p <- ggplot(df_clonal) +
#     geom_point(aes(x=rank, y=n)) +
#     scale_x_continuous(trans="log10") +
#     scale_y_continuous(trans="log10") +
#     theme_bw()
# ggsave("clonal_zipf.png", p)

p <- ggplot(df_clonal) +
    geom_histogram(aes(x=n), bins = 400) +
    # scale_y_continuous(trans="log2") +
    # scale_x_continuous(limits = c(0, 500)) +
    # scale_y_continuous(limits = c(0, 1000)) +
    ylim(0, 500) + xlim(0, 500) +
    theme_bw()
ggsave("clonal_hist.png", p)

df_avj_count <- df_pq %>%
    dplyr::group_by(trav,traj) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = trav, values_from = n) %>%
    as.data.frame()

row.names(df_avj_count) <- df_avj_count$traj
df_avj_count$traj <- NULL

pheatmap(
    log10(df_avj_count),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    filename="avj_count.pdf"
)

df_bvj_count <- df_pq %>%
    dplyr::group_by(trbv,trbj) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = trbv, values_from = n) %>%
    as.data.frame()

row.names(df_bvj_count) <- df_bvj_count$trbj
df_bvj_count$trbj <- NULL

pheatmap(
    log10(df_bvj_count),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    filename="bvj_count.pdf"
)

var_df <- readr::read_tsv("variants.tsv")
var_df_wide <- var_df %>%
    dplyr::select(START, STOP, REF, ALT, CHR, caller) %>%
    dplyr::mutate(present=TRUE) %>%
    tidyr::pivot_wider(names_from = caller, values_from = present)

p <- upset(var_df_wide, unique(var_df$caller), name = "caller")
ggsave("var_upset.pdf", p)