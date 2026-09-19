# Helpers for the integration tests in test_integration_real_data.R.
#
# The integration tests replay the real-data analyses reported in the levi
# supplement. They run only when LEVI_INTEGRATION is set, and the GEO and
# Kang 2018 cases additionally need LEVI_REAL_DATA to point at a directory
# holding the raw files (see inst/scripts/README.md). Nothing here is used
# by the unit tests.

skip_unless_integration <- function() {
    testthat::skip_if_not(nzchar(Sys.getenv("LEVI_INTEGRATION")),
        "set LEVI_INTEGRATION=1 to run the real-data integration tests")
}

real_data_file <- function(...) {
    dir <- Sys.getenv("LEVI_REAL_DATA")
    testthat::skip_if_not(nzchar(dir) && dir.exists(dir),
        "LEVI_REAL_DATA must point at the directory with the raw data files")
    f <- file.path(dir, ...)
    testthat::skip_if_not(file.exists(f), paste("missing", basename(f)))
    f
}

read_saved_network <- function(tag) {
    # STRING or KEGG networks saved as TSV (nodes: name/x/y, edges: V1/V2).
    list(nodes = read.delim(real_data_file(paste0(tag, "_nodes.tsv"))),
         edges = read.delim(real_data_file(paste0(tag, "_edges.tsv"))))
}

# GSE10072 (Landi 2008): HG-U133A RMA series matrix + GPL96 annotation.
# Returns the log2 expression matrix, sample groups and per-gene limma table
# (most significant probe per symbol), exactly as in the supplement scripts.
load_gse10072 <- function() {
    testthat::skip_if_not_installed("limma")
    lines <- readLines(gzfile(real_data_file("GSE10072_series_matrix.txt.gz")))
    hdr <- lines[startsWith(lines, "!")]
    tbl <- read.delim(text = lines[!startsWith(lines, "!") & nzchar(lines)],
        check.names = FALSE, row.names = 1)
    expr <- as.matrix(tbl[rownames(tbl) != "ID_REF", , drop = FALSE])
    field <- function(prefix)
        gsub("\"", "", strsplit(hdr[startsWith(hdr, prefix)][1], "\t")[[1]][-1])
    src <- field("!Sample_source_name_ch1"); ttl <- field("!Sample_title")
    smoke <- sub("Cigarette Smoking Status: ", "",
        gsub("\"", "", strsplit(hdr[grepl("Smoking", hdr)][1], "\t")[[1]][-1]))
    smoke <- c("Current Smoker" = "Current", "Former Smoker" = "Former",
               "Never Smoked" = "Never")[smoke]
    tissue <- ifelse(grepl("Adenocarcinoma|Tumor", src, ignore.case = TRUE),
        "Tumor", "Normal")
    group <- setNames(tissue, colnames(expr))
    patient <- setNames(sub(".*_(GT[0-9]+)\"?$", "\\1", ttl), colnames(expr))
    strata <- setNames(paste(tissue, smoke, sep = "_"), colnames(expr))
    ann <- read.delim(gzfile(real_data_file("GPL96.annot.gz")), skip = 27,
        check.names = FALSE, quote = "")
    ann <- ann[nzchar(ann$`Gene symbol`) & !grepl("///", ann$`Gene symbol`), ]
    sym <- setNames(ann$`Gene symbol`, ann$ID)
    collapse <- function(tt) {
        tt$Symbol <- sym[rownames(tt)]
        tt <- tt[!is.na(tt$Symbol), ]
        tt[!duplicated(tt$Symbol), ]
    }
    design <- model.matrix(~ factor(group, levels = c("Normal", "Tumor")))
    fit <- limma::eBayes(limma::lmFit(expr, design))
    tt <- collapse(limma::topTable(fit, coef = 2, number = Inf, sort.by = "P"))
    list(expr = expr, group = group, patient = patient, strata = strata,
         tt = tt, collapse = collapse)
}

# Per-gene expression matrix for the genes of a network: one probe per gene
# (the most significant one in `tt`).
gene_matrix <- function(expr, tt, genes) {
    probe_of <- setNames(rownames(tt), tt$Symbol)
    gm <- expr[probe_of[genes], , drop = FALSE]
    rownames(gm) <- genes
    gm
}

# Kang 2018 (muscData Kang18_8vs8 as SingleCellExperiment, saved as RDS).
load_kang18 <- function() {
    testthat::skip_if_not_installed("SingleCellExperiment")
    testthat::skip_if_not_installed("edgeR")
    testthat::skip_if_not_installed("limma")
    sce <- readRDS(real_data_file("kang18.rds"))
    sce <- sce[, !is.na(sce$cell)]
    counts <- SingleCellExperiment::counts(sce)
    rownames(counts) <- sub("_ENSG.*$", "", rownames(counts))
    counts <- counts[!duplicated(rownames(counts)), ]
    list(counts = counts, donor = as.character(sce$ind),
         condition = as.character(sce$stim), cell_type = as.character(sce$cell))
}

# Pseudobulk limma-voom table (paired by donor) for one cell type.
kang_pseudobulk_table <- function(k, cell_type) {
    sel <- k$cell_type == cell_type
    pb <- leviPseudobulk(k$counts[, sel], k$donor[sel], k$cell_type[sel],
        k$condition[sel])
    y <- edgeR::DGEList(pb$counts)
    keep <- edgeR::filterByExpr(y, group = pb$condition)
    y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = FALSE])
    design <- model.matrix(~ factor(pb$donor) +
        factor(pb$condition, levels = c("ctrl", "stim")))
    fit <- limma::eBayes(limma::lmFit(limma::voom(y, design), design))
    tt <- limma::topTable(fit, coef = ncol(design), number = Inf, sort.by = "P")
    tt$Symbol <- rownames(tt)
    tt
}

landscape_args <- list(fileTypeInput = "stg", logfc_k = 0.7,
    resolutionValueInput = 30, smoothValueInput = 40, n_perm = 999)

region_p <- function(summary, region) summary$PSpatial[summary$Region == region]
