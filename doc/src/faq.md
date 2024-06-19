# FAQs

## Experimental Design

- **Q: I got scRNA-Seq data in loom/h5ad/others. How can I convert it to a format YASIM-scTCR supports?**

  A: If your format is supported by AnnData, you may convert it using `convert_anndata` command.

- **Q: Do YASIM-scTCR support raw FASTQ reads from manufacturers like 10xGenomics, etc?**

  A: No. The generated FASTQ files should be analyzed using standard pipelines on a cell-to-cell basis.

## Reproducibility

- **Q: We observed that other databases are also providing nucleotide and/or amino-acid sequences for TCR gene segments.**

  A: You are free to use TCR gene segments from [IMGT GENEDB](https://www.imgt.org/download/GENE-DB/), [UniProtKb-SWISSPROT](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz), [IgBLAST Reference TCR](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/). However, since our data is primarily based on CellRanger statistics, we still recommend you to use the same version of TCR segments used by us (i.e., Ensembl release 97). Using a different version of TCR segments will lead to the risk of:

  - Different nomenclatures. The IMGT reference may include `/` in their gene names, which must be removed before further processing.
  - Divergent sequences. The IMGT, SwissProt and Ensembl release 97 have different nucleotide and amino-acid sequences for the same TCR gene segment. IMGT VQUEST may include multiple sequences from different databases for the same gene.

- **Q: What's the version of software used in your studies?**
  
  A: See the following table:

  | Software                                                                                           | Version   |
  |----------------------------------------------------------------------------------------------------|-----------|
  | [GNU Bash](https://www.gnu.org/software/bash)                                                      | 5.2.21(1) |
  | [GNU Grep](https://www.gnu.org/software/grep)                                                      | 3.11      |
  | [GNU Wget](https://www.gnu.org/software/wget)                                                      | 1.21.4    |
  | [GNU Sed](https://www.gnu.org/software/sed)                                                        | 4.9       |
  | [yasim](https://pypi.org/project/yasim/)                                                           | 3.2.1     |
  | [yasim\_sctcr](https://pypi.org/project/yasim-sctcr/)                                              | 1.0.0     |
  | [art\_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) | 2.5.8     |
  | [samtools](https://www.htslib.org)                                                                 | 1.13      |
  | [seqkit](https://bioinf.shenwei.me/seqkit)                                                         | 2.3.0     |
  | [seqtk](https://github.com/lh3/seqtk)                                                              | 1.4-r122  |

- **Q: What's the version of R packages used in your studies?**

  A: See the following code snippet:

    ```r
    library(biomaRt)
    library(tidyverse)
    library(Seurat)
    library(scDesign2)
    library(arrow)
    
    sessionInfo()
    ```

    ```text
    R version 4.3.2 (2023-10-31)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Kali GNU/Linux Rolling
    
    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-openmp/libblas.so.3 
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-openmp/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    
    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
     [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    
    time zone: Asia/Shanghai
    tzcode source: system (glibc)
    
    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] arrow_15.0.1       scDesign2_0.1.0    Seurat_5.0.3       SeuratObject_5.0.1 sp_2.1-3          
     [6] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2       
    [11] readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.0      tidyverse_2.0.0   
    [16] biomaRt_2.58.2    
    
    loaded via a namespace (and not attached):
      [1] RColorBrewer_1.1-3      jsonlite_1.8.8          magrittr_2.0.3          spatstat.utils_3.0-4   
      [5] zlibbioc_1.48.2         vctrs_0.6.5             ROCR_1.0-11             spatstat.explore_3.2-7 
      [9] memoise_2.0.1           RCurl_1.98-1.14         htmltools_0.5.8.1       progress_1.2.3         
     [13] curl_5.2.1              sctransform_0.4.1       parallelly_1.37.1       KernSmooth_2.23-22     
     [17] htmlwidgets_1.6.4       ica_1.0-3               plyr_1.8.9              plotly_4.10.4          
     [21] zoo_1.8-12              cachem_1.0.8            igraph_2.0.3            mime_0.12              
     [25] lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.6-5            R6_2.5.1               
     [29] fastmap_1.1.1           GenomeInfoDbData_1.2.11 fitdistrplus_1.1-11     future_1.33.2          
     [33] shiny_1.8.1.1           digest_0.6.35           colorspace_2.1-0        patchwork_1.2.0        
     [37] AnnotationDbi_1.64.1    S4Vectors_0.40.2        tensor_1.5              RSpectra_0.16-1        
     [41] irlba_2.3.5.1           RSQLite_2.3.6           filelock_1.0.3          progressr_0.14.0       
     [45] spatstat.sparse_3.0-3   fansi_1.0.6             timechange_0.3.0        polyclip_1.10-6        
     [49] abind_1.4-5             httr_1.4.7              compiler_4.3.2          bit64_4.0.5            
     [53] withr_3.0.0             DBI_1.2.2               fastDummies_1.7.3       MASS_7.3-60.0.1        
     [57] rappdirs_0.3.3          tools_4.3.2             lmtest_0.9-40           httpuv_1.6.15          
     [61] future.apply_1.11.2     goftest_1.2-3           glue_1.7.0              nlme_3.1-164           
     [65] promises_1.3.0          grid_4.3.2              Rtsne_0.17              cluster_2.1.6          
     [69] reshape2_1.4.4          generics_0.1.3          spatstat.data_3.0-4     gtable_0.3.4           
     [73] tzdb_0.4.0              data.table_1.15.4       hms_1.1.3               xml2_1.3.6             
     [77] utf8_1.2.4              XVector_0.42.0          spatstat.geom_3.2-9     BiocGenerics_0.48.1    
     [81] RcppAnnoy_0.0.22        ggrepel_0.9.5           RANN_2.6.1              pillar_1.9.0           
     [85] spam_2.10-0             RcppHNSW_0.6.0          later_1.3.2             splines_4.3.2          
     [89] BiocFileCache_2.10.2    lattice_0.22-6          deldir_2.0-4            survival_3.5-8         
     [93] bit_4.0.5               tidyselect_1.2.1        Biostrings_2.70.3       miniUI_0.1.1.1         
     [97] pbapply_1.7-2           gridExtra_2.3           IRanges_2.36.0          scattermore_1.2        
    [101] stats4_4.3.2            Biobase_2.62.0          matrixStats_1.2.0       stringi_1.8.3          
    [105] lazyeval_0.2.2          codetools_0.2-19        cli_3.6.2               uwot_0.1.16            
    [109] xtable_1.8-4            reticulate_1.35.0       munsell_0.5.1           pscl_1.5.9             
    [113] Rcpp_1.0.12             GenomeInfoDb_1.38.8     spatstat.random_3.2-3   globals_0.16.3         
    [117] dbplyr_2.5.0            png_0.1-8               XML_3.99-0.16.1         parallel_4.3.2         
    [121] assertthat_0.2.1        blob_1.2.4              prettyunits_1.2.0       dotCall64_1.1-1        
    [125] bitops_1.0-7            listenv_0.9.1           viridisLite_0.4.2       scales_1.3.0           
    [129] ggridges_0.5.6          leiden_0.4.3.1          crayon_1.5.2            rlang_1.1.3            
    [133] cowplot_1.1.3           KEGGREST_1.42.0   
    ```
  
    Please note that we used [OpenBLAS](https://www.openblas.net/) instead of [Intel oneAPI MKL BLAS](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) since several errors arose while using the latter.
