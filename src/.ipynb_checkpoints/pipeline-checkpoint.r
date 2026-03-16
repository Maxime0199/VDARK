# =============================================================================
# pipeline.R
# Reference-free somatic SNV detection from tumour/normal WGS data.
#
# Workflow:
#   1. Load tumour-specific reads (pre-filtered by extract_tumor_specific_kmers.sh)
#   2. Build read-kmer association table (Python)
#   3. Cluster reads by shared k-mers
#   4. Per cluster: assemble tumour & normal contigs, align, detect SNVs
#   5. Map normal contigs to reference genome and report genomic coordinates
#
# Prerequisites:
#   - extract_tumor_specific_kmers.sh must have been run
#   - Working directory: /srv/home/mlef0011/VDARK/
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(Biostrings)
    library(data.table)
    library(igraph)
    library(pwalign)
    library(Rsamtools)
    library(parallel)
})

source("/srv/home/mlef0011/VDARK/src/utils.r")

WD        <- "/srv/home/mlef0011/VDARK"
MINIMAP   <- "/srv/home/mlef0011/anaconda3/envs/VDARK/bin/minimap2"
REF       <- "/srv/home/mlef0011/rawdata/ref_genome/reference_genome_GRCh37.fa"
N_CORES   <- 20


# =============================================================================
# 1. Load tumour-specific reads
# =============================================================================

log_msg("Loading tumour-specific reads...")

load_reads <- function(r1, r2) {
    fq1 <- readDNAStringSet(r1, format = "fastq")
    fq2 <- readDNAStringSet(r2, format = "fastq")
    data.frame(
        read_ID  = c(names(fq1), names(fq2)),
        sequence = as.character(c(fq1, fq2)),
        stringsAsFactors = FALSE
    )
}

tumour_reads <- load_reads(
    file.path(WD, "rawdata/reads/tumour_R1_tumour_specific.fq"),
    file.path(WD, "rawdata/reads/tumour_R2_tumour_specific.fq")
)
log_msg(nrow(tumour_reads), " tumour-specific reads")


# =============================================================================
# 2. Build read-kmer table
# =============================================================================

log_msg("Building read-kmer table...")
system("/srv/home/mlef0011/anaconda3/condabin/conda run -n VDARK python3 /srv/home/mlef0011/VDARK/src/build_read_kmer_table.py")


# =============================================================================
# 3. Cluster reads by shared k-mers
# =============================================================================

log_msg("Clustering reads...")

read_kmer_df <- fread(file.path(WD, "rawdata/reads/read_kmer_association.tsv"))
read_kmer_df <- read_kmer_df[, if (.N <= 50) .SD, by = kmer]  # drop high-freq k-mers

pairs <- read_kmer_df[, {
    ids <- read_ID
    if (length(ids) >= 2) {
        idx <- combn(length(ids), 2)
        data.table(read_ID.x = ids[idx[1,]], read_ID.y = ids[idx[2,]])
    }
}, by = kmer][, .(weight = .N), by = .(read_ID.x, read_ID.y)]

G        <- graph_from_data_frame(pairs, directed = FALSE)
E(G)$weight <- pairs$weight
clusters <- get("components", envir = asNamespace("igraph"))(G)
log_msg(clusters$no, " clusters")


# =============================================================================
# 4. Per-cluster SNV detection
# =============================================================================

process_cluster <- function(cluster_ID) {

    # Buffer all messages — flush as one block at the end
    logs <- character(0)
    buf  <- function(...) logs <<- c(logs, paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...))

    buf("===== Cluster ", cluster_ID, " =====")

    reads_in  <- names(clusters$membership)[clusters$membership == cluster_ID]
    kmers_sig <- unique(read_kmer_df$kmer[read_kmer_df$read_ID %in% reads_in])
    all_kmers <- unique(unlist(lapply(
        tumour_reads$sequence[tumour_reads$read_ID %in% reads_in], get_kmers
    )))
    buf("  reads: ", length(reads_in), " | specific k-mers: ", length(kmers_sig))

    # ── Assemble tumour contig ─────────────────────────────────────────────────

    contig_tumour <- assemble_kmers(kmers_sig, k = K)
    if (is.null(contig_tumour) || nchar(contig_tumour[1]) < K) {
        buf("  Tumour contig too short -- skipping")
        message(paste(logs, collapse = "\n")); return(NULL)
    }
    buf("  Tumour contig: ", nchar(contig_tumour[1]), " bp")

    # ── Fetch normal reads via flanking k-mers ────────────────────────────────

    flanking_kmers <- setdiff(all_kmers, kmers_sig)
    flanking_kmers <- flanking_kmers[grepl("^[ACGT]+$", flanking_kmers)]
    tmp_cl      <- file.path(WD, "tmp", cluster_ID)
    dir.create(tmp_cl, showWarnings = FALSE, recursive = TRUE)
    flanking_fa <- file.path(tmp_cl, paste0("flanking_", cluster_ID, ".fa"))
    kmc_db      <- file.path(tmp_cl, paste0("flanking_kmc_", cluster_ID))

    writeLines(paste0(">k_", seq_along(flanking_kmers), "\n", flanking_kmers), flanking_fa)

    system(paste(KMC, "-k31 -t12 -ci1 -fm", flanking_fa, kmc_db, tmp_cl),
           ignore.stdout = TRUE, ignore.stderr = TRUE)

    NORMAL_R1 <- file.path(WD, "rawdata/reads/normal_R1.fq")
    NORMAL_R2 <- file.path(WD, "rawdata/reads/normal_R2.fq")
    fq_R1_out <- file.path(tmp_cl, paste0("normal_R1_locus_", cluster_ID, ".fq"))
    fq_R2_out <- file.path(tmp_cl, paste0("normal_R2_locus_", cluster_ID, ".fq"))

    system(paste(KMC_TOOLS, "filter", kmc_db, "-ci1", NORMAL_R1, fq_R1_out),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste(KMC_TOOLS, "filter", kmc_db, "-ci1", NORMAL_R2, fq_R2_out),
           ignore.stdout = TRUE, ignore.stderr = TRUE)

    normal_df <- tryCatch(
        load_reads(fq_R1_out, fq_R2_out),
        error = function(e) data.frame(read_ID = character(), sequence = character())
    )

    # Precise filtering: keep reads with >= 30 flanking k-mer hits
    if (nrow(normal_df) > 0) {
        pdict  <- PDict(unique(c(flanking_kmers,
                                 as.character(reverseComplement(DNAStringSet(flanking_kmers))))))
        n_hits <- vapply(vwhichPDict(pdict, DNAStringSet(normal_df$sequence)),
                         length, integer(1))
        normal_df <- normal_df[n_hits >= 30, ]
    }

    if (nrow(normal_df) == 0) {
        buf("  No normal reads at locus -- skipping")
        message(paste(logs, collapse = "\n")); return(NULL)
    }
    buf("  Normal reads retained: ", nrow(normal_df))

    # ── Assemble normal contig ─────────────────────────────────────────────────

    kmer_counts      <- table(unlist(lapply(normal_df$sequence, get_kmers)))
    buf("  Normal k-mers with cov>=10: ", sum(kmer_counts >= 10), " / ", length(kmer_counts))
    kmer_to_assemble <- names(kmer_counts[kmer_counts >= 10])

    if (length(kmer_to_assemble) == 0) {
        buf("  No k-mer with cov>=10 in normal -- skipping")
        message(paste(logs, collapse = "\n")); return(NULL)
    }

    contig_normal <- assemble_kmers(kmer_to_assemble, k = K)
    if (is.null(contig_normal) || nchar(contig_normal[1]) < K) {
        buf("  Normal contig too short -- skipping")
        message(paste(logs, collapse = "\n")); return(NULL)
    }
    buf("  Normal contig: ", nchar(contig_normal[1]), " bp")

    # ── Align tumour vs normal ─────────────────────────────────────────────────

    aln <- best_alignment(contig_normal, contig_tumour)
    format_alignment(aln, cluster_ID, buf = buf)

    mm <- mismatchTable(aln)
    if (nrow(mm) == 0) {
        buf("  No mismatches detected")
        message(paste(logs, collapse = "\n")); return(NULL)
    }

    # ── Germline filter ────────────────────────────────────────────────────────

    germline   <- is_germline(mm, contig_normal, cluster_ID, buf = buf)
    mm_somatic <- mm[!germline, ]

    if (nrow(mm_somatic) == 0) {
        buf("  No somatic variants (all germline)")
        message(paste(logs, collapse = "\n")); return(NULL)
    }
    buf("  Somatic SNV(s): ", nrow(mm_somatic))

    # Flush all buffered messages as one block
    message(paste(logs, collapse = "\n"))
    unlink(tmp_cl, recursive = TRUE)
    list(
        cluster_ID    = cluster_ID,
        mismatches    = mm_somatic,
        contig_tumour = contig_tumour,
        contig_normal = contig_normal
    )
}

SNV_list <- Filter(Negate(is.null),
                   mclapply(seq_len(clusters$no), process_cluster, mc.cores = N_CORES))

log_msg(length(SNV_list), " cluster(s) with somatic SNVs")


# =============================================================================
# 5. Map normal contigs to reference & report genomic coordinates
# =============================================================================

log_msg("Mapping contigs to reference genome...")

contig_fa  <- file.path(WD, "rawdata/contigs/contigs_normal.fa")
contig_sam <- file.path(WD, "rawdata/contigs/contigs_normal_aln.sam")

writeLines(
    unlist(lapply(SNV_list, function(x)
        c(paste0(">contig_cluster_", x$cluster_ID), x$contig_normal[1])
    )),
    contig_fa
)

system(paste(MINIMAP, "-ax sr", REF, contig_fa, ">", contig_sam))

sam <- read.table(contig_sam, sep = "\t", comment.char = "@", fill = TRUE)
# Filter SAM to keep only one line per contig (first occurrence of each QNAME)
sam_clean <- sam[!duplicated(sam[, 1]), ]
# Then reorder to match SNV_list order
contig_names <- paste0("contig_cluster_", sapply(SNV_list, function(x) x$cluster_ID))
sam_matched  <- sam_clean[match(contig_names, sam_clean[, 1]), ]
# ── Build SNV table ───────────────────────────────────────────────────────────

snv_rows <- lapply(seq_along(SNV_list), function(i) {
    x    <- SNV_list[[i]]
    mm   <- x$mismatches
    flag <- as.integer(sam_matched[i, 2])

    if (is.na(flag) || bitwAnd(flag, 4) != 0 || bitwAnd(flag, 2048) != 0) {
        log_msg("  Cluster ", x$cluster_ID, " | contig unmapped or supplementary -- skipping")
        return(NULL)
    }

    is_reverse <- bitwAnd(flag, 16) != 0
    chrom      <- as.character(sam_matched[i, 3])
    pos        <- as.integer(sam_matched[i, 4])
    contig_len <- nchar(as.character(sam_matched[i, 10]))

    cigar     <- as.character(sam_matched[i, 6])
    m         <- regmatches(cigar, regexpr("^([0-9]+)[SH]", cigar))
    left_clip <- if (length(m) > 0) as.integer(sub("[SH]", "", m)) else 0

    lapply(seq_len(nrow(mm)), function(j) {
        mut_pos <- if (is_reverse) {
            pos + contig_len - mm$PatternStart[j]
        } else {
            pos + mm$PatternStart[j] - 1L
        }
        ref_nt <- as.character(mm$PatternSubstring[j])
        alt_nt <- as.character(mm$SubjectSubstring[j])
        if (is_reverse) {
            ref_nt <- as.character(reverseComplement(DNAString(ref_nt)))
            alt_nt <- as.character(reverseComplement(DNAString(alt_nt)))
        }
        data.frame(CHROM = chrom, POS = mut_pos, REF = ref_nt, ALT = alt_nt,
                   stringsAsFactors = FALSE)
    })
})

snv_df <- do.call(rbind, Filter(Negate(is.null), unlist(snv_rows, recursive = FALSE)))
snv_df <- do.call(rbind, unlist(snv_rows, recursive = FALSE))
snv_df <- snv_df[order(snv_df$CHROM, snv_df$POS), ]

# ── Print ─────────────────────────────────────────────────────────────────────

cat("\n=================== Somatic SNVs ===================\n")
cat("#CHROM\tPOS\tREF\tALT\n")
for (i in seq_len(nrow(snv_df))) {
    cat(sprintf("%s\t%d\t%s\t%s\n",
                snv_df$CHROM[i], snv_df$POS[i], snv_df$REF[i], snv_df$ALT[i]))
}

# ── Write TSV ─────────────────────────────────────────────────────────────────

tsv_out <- file.path(WD, "results/somatic_snvs.tsv")
dir.create(dirname(tsv_out), showWarnings = FALSE)
writeLines(paste(c("#CHROM", "POS", "REF", "ALT"), collapse = "\t"), tsv_out)
write.table(snv_df, tsv_out, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE, append = TRUE)
log_msg(nrow(snv_df), " SNV(s) written to ", tsv_out)