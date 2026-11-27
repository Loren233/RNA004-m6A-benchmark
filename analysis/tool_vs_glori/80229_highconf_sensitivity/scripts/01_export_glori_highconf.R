#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(GenomicRanges)
})

input_rds <- file.path(Sys.getenv("RNA004_STAT_ROOT", unset = "analysis"), "ground_truth/m6aconquer/raw/GLORI_hg38_WT_m6A_SE.rds")

output_tsv <- file.path(Sys.getenv("RNA004_STAT_ROOT", unset = "analysis"), "tool_vs_glori/80229_highconf_sensitivity/results/GLORI_HEK_HeLa.highconf.tsv.gz")

se <- readRDS(input_rds)

cat("Loaded:", input_rds, "\n")
cat("Dimensions:", nrow(se), "x", ncol(se), "\n")

# ------------------------------------------------------------
# Find the four samples used in our main GLORI comparison
# ------------------------------------------------------------

find_sample <- function(se, accession) {

    cn <- colnames(se)

    hit <- which(grepl(accession, cn, fixed = TRUE))

    if (length(hit) == 1) {
        return(hit)
    }

    cd <- as.data.frame(colData(se))

    for (nm in colnames(cd)) {
        vals <- as.character(cd[[nm]])
        hit <- which(grepl(accession, vals, fixed = TRUE))

        if (length(hit) == 1) {
            return(hit)
        }
    }

    stop(
        paste(
            "Could not uniquely locate sample:",
            accession
        )
    )
}

# Main GT selected previously
hek_rep1  <- find_sample(se, "SRR21356250")
hek_rep2  <- find_sample(se, "SRR21356251")

hela_rep1 <- find_sample(se, "SRR21356245")
hela_rep2 <- find_sample(se, "SRR21356246")

cat("HEK rep1 :", colnames(se)[hek_rep1], "\n")
cat("HEK rep2 :", colnames(se)[hek_rep2], "\n")
cat("HeLa rep1:", colnames(se)[hela_rep1], "\n")
cat("HeLa rep2:", colnames(se)[hela_rep2], "\n")

# ------------------------------------------------------------
# Extract assays
# ------------------------------------------------------------

m6a   <- assay(se, "m6A")
total <- assay(se, "Total")
padj  <- assay(se, "AdjPvalue")

# ------------------------------------------------------------
# Pool biological replicates by counts
# ------------------------------------------------------------

hek_m6a <-
    as.numeric(m6a[, hek_rep1]) +
    as.numeric(m6a[, hek_rep2])

hek_total <-
    as.numeric(total[, hek_rep1]) +
    as.numeric(total[, hek_rep2])

hela_m6a <-
    as.numeric(m6a[, hela_rep1]) +
    as.numeric(m6a[, hela_rep2])

hela_total <-
    as.numeric(total[, hela_rep1]) +
    as.numeric(total[, hela_rep2])

hek_beta <- ifelse(
    hek_total > 0,
    hek_m6a / hek_total,
    NA_real_
)

hela_beta <- ifelse(
    hela_total > 0,
    hela_m6a / hela_total,
    NA_real_
)

delta_beta <- hek_beta - hela_beta

# ------------------------------------------------------------
# Replicate-level adjusted p-values
#
# min = at least one replicate significant
# max = both replicates significant
# ------------------------------------------------------------

pmin_na <- function(a, b) {

    out <- pmin(a, b, na.rm = TRUE)

    both_na <- is.na(a) & is.na(b)
    out[both_na] <- NA_real_

    out
}

pmax_na <- function(a, b) {

    out <- pmax(a, b, na.rm = TRUE)

    both_na <- is.na(a) & is.na(b)
    out[both_na] <- NA_real_

    out
}

hek_p1 <- as.numeric(padj[, hek_rep1])
hek_p2 <- as.numeric(padj[, hek_rep2])

hela_p1 <- as.numeric(padj[, hela_rep1])
hela_p2 <- as.numeric(padj[, hela_rep2])

hek_p_min  <- pmin_na(hek_p1, hek_p2)
hek_p_max  <- pmax_na(hek_p1, hek_p2)

hela_p_min <- pmin_na(hela_p1, hela_p2)
hela_p_max <- pmax_na(hela_p1, hela_p2)

# ------------------------------------------------------------
# Coordinates
# ------------------------------------------------------------

rr <- rowRanges(se)

out <- data.frame(
    chr = as.character(seqnames(rr)),
    pos_1based = start(rr),
    strand = as.character(strand(rr)),

    glori_m6a_hek = hek_m6a,
    glori_total_hek = hek_total,
    glori_beta_hek = hek_beta,

    glori_m6a_hela = hela_m6a,
    glori_total_hela = hela_total,
    glori_beta_hela = hela_beta,

    glori_delta_beta = delta_beta,

    glori_adj_p_hek_rep1 = hek_p1,
    glori_adj_p_hek_rep2 = hek_p2,
    glori_adj_p_hek_min = hek_p_min,
    glori_adj_p_hek_max = hek_p_max,

    glori_adj_p_hela_rep1 = hela_p1,
    glori_adj_p_hela_rep2 = hela_p2,
    glori_adj_p_hela_min = hela_p_min,
    glori_adj_p_hela_max = hela_p_max,

    stringsAsFactors = FALSE
)

cat("Rows:", nrow(out), "\n")

# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

if (requireNamespace("data.table", quietly = TRUE)) {

    data.table::fwrite(
        out,
        file = output_tsv,
        sep = "\t"
    )

} else {

    con <- gzfile(output_tsv, "wt")

    write.table(
        out,
        file = con,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    close(con)
}

cat("Written:", output_tsv, "\n")
