# =============================================================================
# filter_tumour_specific_kmers.r
# Combine Tier 1 (normal=0, auto-keep) et Tier 2 (zone grise, test binomial)
# en une seule liste finale de k-mers tumeur-spécifiques.
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
    make_option("--histo_tumour",  type = "character"),
    make_option("--histo_normal",  type = "character"),
    make_option("--tier1",         type = "character"),
    make_option("--grey_tumour",   type = "character"),
    make_option("--grey_normal",   type = "character"),
    make_option("--p_threshold",   type = "double", default = 0.01),
    make_option("--out",           type = "character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

log_msg <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

# ── Profondeurs globales -> p_null ────────────────────────────────────────────
histo_tumour <- read.table(opt$histo_tumour, col.names = c("count", "n_kmers"))
histo_normal <- read.table(opt$histo_normal, col.names = c("count", "n_kmers"))

N_T <- sum(histo_tumour$count * histo_tumour$n_kmers)
N_N <- sum(histo_normal$count * histo_normal$n_kmers)
p_null <- N_T / (N_T + N_N)

log_msg("N_T=", N_T, " N_N=", N_N, " p_null=", round(p_null, 4))

# ── Tier 2 : test binomial sur la zone grise ─────────────────────────────────
tumour_grey_kmer <- read.table(opt$grey_tumour, col.names = c("kmer", "tumour"))
normal_grey_kmer <- read.table(opt$grey_normal, col.names = c("kmer", "normal"))
grey_kmer <- merge(tumour_grey_kmer, normal_grey_kmer)

grey_kmer$n_total <- grey_kmer$tumour + grey_kmer$normal
grey_kmer$p_value <- pbinom(grey_kmer$tumour - 1, grey_kmer$n_total, p_null, lower.tail = FALSE)
grey_kmer$p_adj   <- p.adjust(grey_kmer$p_value, method = "BH")

tumour_specific_grey <- grey_kmer[grey_kmer$p_adj < opt$p_threshold, c("kmer", "tumour")]
log_msg(nrow(tumour_specific_grey), " / ", nrow(grey_kmer),
        " k-mer(s) significatifs dans la zone grise (p_adj < ", opt$p_threshold, ")")

# ── Tier 1 : garde direct ─────────────────────────────────────────────────────
tier1 <- read.table(opt$tier1, col.names = c("kmer", "tumour"))
log_msg(nrow(tier1), " k-mer(s) en Tier 1 (normal=0)")

# ── Combine et écrit le résultat final ────────────────────────────────────────
tumour_specific_final <- rbind(tier1, tumour_specific_grey)
log_msg(nrow(tumour_specific_final), " k-mer(s) tumeur-spécifiques au total -> ", opt$out)

writeLines(tumour_specific_final$kmer, opt$out)