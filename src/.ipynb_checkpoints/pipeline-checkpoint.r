# =============================================================================
# pipeline.r
# Reference-free somatic SNV detection from tumour/normal WGS data,
# using a De Bruijn Graph (DBG) assembler.
#
# Prerequisites:
#   - extract_tumor_specific_kmers.sh must have been run
#   - Working directory: /srv/home/mlef0011/VDARK/
#   - utils.r and dbg_assembly_igraph.r must live in /srv/home/mlef0011/VDARK/src/
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
source("/srv/home/mlef0011/VDARK/src/dbg_assembly_igraph.r")

WD        <- "/srv/home/mlef0011/VDARK/"
MINIMAP   <- "/srv/home/mlef0011/anaconda3/envs/VDARK/bin/minimap2"
REF       <- "/srv/home/mlef0011/rawdata/ref_genome/reference_genome_GRCh37.fa"
N_CORES   <- 20
K         <- 31
COV       <- 55

KMER_MAX_FREQ          <- 50  # drop high-freq k-mers in tumour specific k-mers before clustering
MIN_SHARED_KMERS       <- 2   # minimum shared k-mers between 2 reads to build a graph edge
NORMAL_KMER_MIN_COV    <- 10
CLUSTER_SIZE_CAP       <- 100
MIN_MEAN_WEIGHT        <- 5
MIN_ALN_SCORE          <- 30
MAX_N_SIG_KMERS_NORMAL <- 5
MAX_N_KMER_SIG         <- K
MIN_N_KMER_SIG        <- 10
GERMLINE_TSH          <- 10
P_NULL_SEQ_ERROR      <- 0.005   # taux d'erreur de sequencage attendu (bruit de fond)
VAF_PVAL_THRESHOLD    <- 0.001   # seuil de significativite classique


label <- paste0("chr21_k", K,"_filtered")

label_file <- sprintf(
  "chr21_k%d_kmf%d_msk%d_ncov%d_cap%d_mw%d_mas%d_maxnsn%d_maxnsk%d_minnsk%d_GERMF%d_filtered",
  K,
  KMER_MAX_FREQ,
  MIN_SHARED_KMERS,
  NORMAL_KMER_MIN_COV,
  CLUSTER_SIZE_CAP,
  MIN_MEAN_WEIGHT,
  MIN_ALN_SCORE,
  MAX_N_SIG_KMERS_NORMAL,
  MAX_N_KMER_SIG,
  MIN_N_KMER_SIG,
  GERMLINE_TSH
)

# ── Per-cluster error log ─────────────────────────────────────────────────────
# Any error inside process_cluster() is caught, appended here, and the
# cluster is skipped (returns NULL) so the rest of the run isn't affected.
ERROR_LOG <- file.path(WD, "logs", paste0("errors_", label_file, ".log"))
dir.create(dirname(ERROR_LOG), showWarnings = FALSE, recursive = TRUE)
if (file.exists(ERROR_LOG)) file.remove(ERROR_LOG)

log_error <- function(cluster_ID, e) {
    cat(sprintf("[%s] Cluster %s ERROR: %s\n",
                format(Sys.time(), "%H:%M:%S"), cluster_ID, conditionMessage(e)),
        file = ERROR_LOG, append = TRUE)
}

# ── Per-cluster stage log ─────────────────────────────────────────────────────
# Every exit point of process_cluster() (early return or success) writes one
# line here: cluster_ID <TAB> stage. Read back at the end of the run to get
# a breakdown of where clusters are lost — much more actionable than a bare
# success/failure count.
STAGE_LOG <- file.path(WD, "logs", paste0("stages_", label_file, ".log"))
dir.create(dirname(STAGE_LOG), showWarnings = FALSE, recursive = TRUE)
if (file.exists(STAGE_LOG)) file.remove(STAGE_LOG)

log_stage <- function(cluster_ID, stage) {
    cat(sprintf("%s\t%s\n", cluster_ID, stage), file = STAGE_LOG, append = TRUE)
}

quiet <- function(...) invisible(NULL)  # used to silence chatty helper functions


# =============================================================================
# Assembly wrapper
# =============================================================================

assembly <- function(kmers, k, kmer_freq, coverage = NULL, mode = c("normal", "tumour")) {
    mode <- match.arg(mode)
    assemble_dbg(kmers, k = k, kmer_freq = kmer_freq, coverage = coverage, mode = mode)
}


# =============================================================================
# 1. Load tumour-specific reads
# =============================================================================

log_msg("Loading tumour-specific reads...")

load_reads <- function(r1, r2) {
    fq1 <- readDNAStringSet(r1, format = "fastq")
    fq2 <- readDNAStringSet(r2, format = "fastq")
    data.frame(
        read_ID  = c(paste0(names(fq1), "/R1"), paste0(names(fq2), "/R2")),
        sequence = as.character(c(fq1, fq2)),
        stringsAsFactors = FALSE
    )
}

tumour_reads <- load_reads(
    file.path(WD, paste0("rawdata/reads/tumour_R1_tumour_specific_", label, ".fq")),
    file.path(WD, paste0("rawdata/reads/tumour_R2_tumour_specific_", label, ".fq"))
)
setDT(tumour_reads)
setkey(tumour_reads, read_ID)
log_msg(nrow(tumour_reads), " tumour-specific reads")


# =============================================================================
# 2. Build read-kmer table
# =============================================================================

log_msg("Building read-kmer table...")
system(paste0("/srv/home/mlef0011/anaconda3/condabin/conda run -n VDARK python3 ",
              WD, "src/build_read_kmer_table.py -k 31 -f ",
              WD, "rawdata/kmer/tumour_specific_", label, ".txt -o ",
              WD, "rawdata/reads/read_kmer_association_", label, ".tsv"))


# =============================================================================
# 3. Cluster reads by shared k-mers
# =============================================================================

log_msg("Clustering reads...")

read_kmer_df <- fread(file.path(WD, paste0("rawdata/reads/read_kmer_association_", label, ".tsv")))
read_kmer_df <- read_kmer_df[, if (.N <= KMER_MAX_FREQ) .SD, by = kmer]

pairs <- read_kmer_df[, {
    ids <- read_ID
    if (length(ids) >= 2) {
        idx <- combn(length(ids), 2)
        data.table(read_ID.x = ids[idx[1, ]], read_ID.y = ids[idx[2, ]])
    }
}, by = kmer][, .(weight = .N), by = .(read_ID.x, read_ID.y)]

pairs <- pairs[weight >= MIN_SHARED_KMERS]

G           <- graph_from_data_frame(pairs, directed = FALSE)
E(G)$weight <- pairs$weight
clusters    <- get("components", envir = asNamespace("igraph"))(G)

log_msg(clusters$no, " clusters")

big_clusters <- which(clusters$csize > CLUSTER_SIZE_CAP)
log_msg(length(big_clusters), " cluster(s) excluded (>", CLUSTER_SIZE_CAP, " reads) — likely repetitive")

cluster_of <- clusters$membership

mean_weight_dt <- pairs[, cluster_ID := cluster_of[read_ID.x]
                        ][, .(mean_weight = mean(weight), n_edges = .N),
                          by = cluster_ID]
setkey(mean_weight_dt, cluster_ID)

get_mean_weight <- function(cid) {
    mw <- mean_weight_dt[.(cid), mean_weight]
    if (length(mw) == 0 || is.na(mw)) return(NA_real_)
    mw
}

# ── Filtre mean_weight, appliqué une fois pour tout le run ───────────────────
all_ids             <- seq_len(clusters$no)
mw_all              <- mean_weight_dt[.(all_ids), mean_weight]  # NA = pas dans pairs (cluster isolé/faible signal)
low_weight_clusters <- all_ids[is.na(mw_all) | mw_all < MIN_MEAN_WEIGHT]

log_msg(length(low_weight_clusters), " cluster(s) excluded (mean_weight < ", MIN_MEAN_WEIGHT, ")")

valid_clusters <- setdiff(all_ids, union(big_clusters, low_weight_clusters))
log_msg(length(valid_clusters), " valid cluster(s) after all pre-filters")

# ── Pre-split read_kmer_df by cluster, once ──────────────────────────────────
read_kmer_df[, cluster_ID := cluster_of[read_ID]]
kmer_by_cluster <- split(read_kmer_df[!is.na(cluster_ID)], by = "cluster_ID")

log_msg("Errors (if any) will be logged to ", ERROR_LOG)
log_msg("Per-cluster stage outcomes will be logged to ", STAGE_LOG)


# =============================================================================
# 4. Per-cluster SNV detection
# =============================================================================

process_cluster <- function(cluster_ID, k = K, idx = NA, n_total = NA) {

    # ── Heartbeat: printed directly to console, ~50 lines over the whole run ──
    if (!is.na(idx) && idx %% PROGRESS_EVERY == 0) {
        elapsed <- round(as.numeric(Sys.time() - PIPELINE_START, units = "secs"))
        cat(sprintf("[%s] %d / %d clusters lances (%ds ecoules)\n",
                    format(Sys.time(), "%H:%M:%S"), idx, n_total, elapsed))
        flush(stdout())
    }

    tmp_cl <- file.path(WD, "tmp", label_file, cluster_ID)
    unlink(tmp_cl, recursive = TRUE)
    dir.create(tmp_cl, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(tmp_cl, recursive = TRUE), add = TRUE)  # always clean up

    tryCatch({

        # mean_weight is already >= MIN_MEAN_WEIGHT since valid_clusters was
        # pre-filtered before mclapply — just retrieve it for the output row.
        mean_weight <- get_mean_weight(cluster_ID)

        reads_in <- names(clusters$membership)[clusters$membership == cluster_ID]

        cl_df <- kmer_by_cluster[[as.character(cluster_ID)]]
        if (is.null(cl_df) || nrow(cl_df) == 0) {
            log_stage(cluster_ID, "no_kmer_data"); return(NULL)
        }
        kmers_sig <- unique(cl_df$kmer)
        
        if (length(kmers_sig) <= MIN_N_KMER_SIG) {
            log_stage(cluster_ID, "too_few_kmer_sig"); return(NULL)
        }

        seqs      <- tumour_reads[reads_in, sequence, nomatch = NULL]
        all_kmers <- unique(unlist(lapply(seqs, get_kmers, k = k)))

        # ── Assemble tumour contig ───────────────────────────────────────────
        sig_counts <- cl_df[, .N, by = kmer]
        kmer_freq  <- setNames(sig_counts$N, sig_counts$kmer)
        asm_t <- assembly(kmers = kmers_sig, k = k, kmer_freq = kmer_freq,
                          coverage = COV, mode = "tumour")
        contig_tumour <- asm_t$contigs
        if (is.null(contig_tumour) || nchar(contig_tumour[1]) < k) {
            log_stage(cluster_ID, "tumour_contig_too_short"); return(NULL)
        }

        # ── Fetch normal reads via flanking k-mers ───────────────────────────
        flanking_kmers <- setdiff(all_kmers, kmers_sig)
        flanking_kmers <- flanking_kmers[grepl("^[ACGT]+$", flanking_kmers)]

        flanking_fa <- file.path(tmp_cl, paste0("flanking_", cluster_ID, "_", label, ".fa"))
        kmc_db      <- file.path(tmp_cl, paste0("flanking_kmc_", label, cluster_ID))

        writeLines(paste0(">k_", seq_along(flanking_kmers), "\n", flanking_kmers), flanking_fa)

        system(paste(KMC, paste0("-k", k), "-t10 -ci1 -fm",
                     flanking_fa, kmc_db, tmp_cl),
               ignore.stdout = TRUE, ignore.stderr = TRUE)

        fq_R1_out <- file.path(tmp_cl, paste0("normal_R1_locus_", cluster_ID, "_", label, ".fq"))
        fq_R2_out <- file.path(tmp_cl, paste0("normal_R2_locus_", cluster_ID, "_", label, ".fq"))

        NORMAL_R1 <- file.path(WD, "rawdata/reads/normal_chr21_R1.fq")
        NORMAL_R2 <- file.path(WD, "rawdata/reads/normal_chr21_R2.fq")

        tsh <- max(5, min(50, round(min(120, length(flanking_kmers)) * 0.25)))

        system(paste(KMC_TOOLS, "-t10", "filter", kmc_db, "-ci1",
                     NORMAL_R1, paste0("-ci", tsh), fq_R1_out),
               ignore.stdout = TRUE, ignore.stderr = TRUE)
        system(paste(KMC_TOOLS, "-t10", "filter", kmc_db, "-ci1",
                     NORMAL_R2, paste0("-ci", tsh), fq_R2_out),
               ignore.stdout = TRUE, ignore.stderr = TRUE)

        normal_df <- tryCatch(
            load_reads(fq_R1_out, fq_R2_out),
            error = function(e) data.frame(read_ID = character(), sequence = character())
        )
        if (nrow(normal_df) == 0) {
            log_stage(cluster_ID, "no_normal_reads"); return(NULL)
        }

        tumour_kmers_rc <- get_kmers(contig_tumour[2], k)  # computed once, not per normal read
        normal_df$n_sig_kmers <- vapply(
            normal_df$sequence,
            function(seq) sum(get_kmers(seq, k) %in% tumour_kmers_rc),
            integer(1)
        )
        n_sig_kmers_in_normal <- sum(normal_df$n_sig_kmers)

        if (n_sig_kmers_in_normal >= MAX_N_SIG_KMERS_NORMAL) {
            log_stage(cluster_ID, "sig_kmers_leak_normal"); return(NULL)
        }

        # ── Assemble normal contig ───────────────────────────────────────────
        kmer_counts <- table(unlist(lapply(normal_df$sequence, get_kmers)))
        n_normal_kmers_total  <- length(kmer_counts)
        n_normal_kmers_cov_ok <- sum(kmer_counts >= NORMAL_KMER_MIN_COV)
        kmer_to_assemble <- names(kmer_counts[kmer_counts >= NORMAL_KMER_MIN_COV])
        if (length(kmer_to_assemble) == 0) {
            log_stage(cluster_ID, "no_normal_kmer_cov"); return(NULL)
        }

        kmer_freq_vec <- kmer_counts[kmer_counts >= NORMAL_KMER_MIN_COV]
        asm_n <- assembly(kmers = kmer_to_assemble, k = k, kmer_freq = kmer_freq_vec,
                          coverage = COV, mode = "normal")
        contig_normal     <- asm_n$contigs
        normal_bubbles_df <- asm_n$bubbles
        if (is.null(contig_normal) || nchar(contig_normal[1]) < k) {
            log_stage(cluster_ID, "normal_contig_too_short"); return(NULL)
        }

        # ── Align tumour vs normal ────────────────────────────────────────────
        aln_result <- best_alignment(contig_normal, contig_tumour)
        aln        <- aln_result$alignment
        normal_rc  <- aln_result$normal_rc
        aln_score  <- score(aln)

        if (aln_score < MIN_ALN_SCORE) {
            log_stage(cluster_ID, "aln_score_too_low"); return(NULL)
        }

        mm <- mismatchTable(aln)
        n_mismatches_total <- nrow(mm)
        if (n_mismatches_total == 0) {
            log_stage(cluster_ID, "no_mismatches"); return(NULL)
        }

        contig_strand <- if (normal_rc) contig_normal[2] else contig_normal[1]

        # ── Germline filter: bubbles (free, DBG) then k-mer presence in normal ─
        bubble_germline <- filter_germline_bubbles(
            mm, contig_strand, normal_bubbles_df, normal_rc, buf = quiet
        )
        mm_for_kmc <- mm[!bubble_germline, ]

        germline_kmc <- is_germline(mm_for_kmc, contig_strand,
                                    normal_df$sequence,
                                    k = k, threshold = GERMLINE_TSH, buf = quiet)

        germline                   <- rep(TRUE, nrow(mm))
        germline[!bubble_germline] <- germline_kmc
        mm_somatic                 <- mm[!germline, ]
        if (nrow(mm_somatic) == 0) {
            log_stage(cluster_ID, "all_germline"); return(NULL)
        }

        
        if (length(kmers_sig) > MAX_N_KMER_SIG & nrow(mm_somatic) ==1) {
            log_stage(cluster_ID, "too_many_kmer_sig"); return(NULL)
        }

        # ── VAF via k-mers + test de significativite binomial ─────────────────
        all_tumour_at_locus <- extract_tumour_reads_at_locus(kmc_db, tmp_cl)
        vaf_df <- estimate_vaf(mm_somatic, contig_strand, all_tumour_at_locus, k = k)

        depth <- vaf_df$n_alt + vaf_df$n_ref
        vaf_df$p_value <- ifelse(
            is.na(vaf_df$n_alt) | is.na(depth) | depth == 0,
            NA_real_,
            pbinom(vaf_df$n_alt - 1, depth, P_NULL_SEQ_ERROR, lower.tail = FALSE)
        )
        vaf_df$is_true_mutation <- !is.na(vaf_df$p_value) & vaf_df$p_value < VAF_PVAL_THRESHOLD

        keep       <- vaf_df$is_true_mutation
        mm_somatic <- mm_somatic[keep, , drop = FALSE]
        vaf_df     <- vaf_df[keep, , drop = FALSE]

        if (nrow(mm_somatic) == 0) {
            log_stage(cluster_ID, "no_true_mutation"); return(NULL)
        }

        # ── Alignment quality metrics ─────────────────────────────────────────
        len_tumour <- nchar(contig_tumour[1])
        len_normal <- nchar(contig_normal[1])
        aln_score_norm <- aln_score / min(len_tumour, len_normal)
        pid_aln <- tryCatch(pid(aln, type = "PID1"), error = function(e) NA_real_)

        log_stage(cluster_ID, "OK")

        list(
            cluster_ID    = cluster_ID,
            mismatches    = mm_somatic,
            vaf_df        = vaf_df,
            contig_tumour = contig_tumour,
            length_contig_tumour = len_tumour,
            contig_normal = contig_normal,
            length_contig_normal = len_normal,
            normal_rc     = normal_rc,
            n_kmer_sig    = length(kmers_sig),
            n_kmer_sig_dist_k = abs(length(kmers_sig) - K),
            n_reads       = length(reads_in),
            n_normal_reads = nrow(normal_df),
            n_normal_kmers_total  = n_normal_kmers_total,
            n_normal_kmers_cov_ok = n_normal_kmers_cov_ok,
            frac_normal_kmers_cov_ok = n_normal_kmers_cov_ok / n_normal_kmers_total,
            n_mismatches_total   = n_mismatches_total,
            n_mismatches_somatic = nrow(mm_somatic),
            mean_weight   = mean_weight,
            aln_score     = aln_score,
            aln_score_norm = aln_score_norm,
            pid_aln       = pid_aln,
            n_sig_kmers_in_normal = n_sig_kmers_in_normal
        )

    }, error = function(e) {
        log_error(cluster_ID, e)
        log_stage(cluster_ID, "ERROR")
        NULL
    })
}

setDTthreads(1)
Sys.setenv(OMP_NUM_THREADS = "1")

#valid_clusters = valid_clusters[1:100]
             
PIPELINE_START <- Sys.time()
N_VALID        <- length(valid_clusters)
PROGRESS_EVERY <- max(1L, round(N_VALID / 50))   # ~50 heartbeat lines over the whole run

log_msg("Lancement sur ", N_VALID, " clusters, heartbeat tous les ", PROGRESS_EVERY, " clusters")

system.time(raw_results <- mclapply(seq_along(valid_clusters), function(i) {
    process_cluster(valid_clusters[i], idx = i, n_total = N_VALID)
}, mc.cores = N_CORES))

SNV_list <- Filter(function(x) is.list(x) && !is.null(x$cluster_ID), raw_results)

n_errors <- sum(vapply(raw_results, function(x) inherits(x, "try-error"), logical(1)))
if (n_errors > 0)
    log_msg("WARNING: ", n_errors, " cluster(s) failed at the mclapply level — see ", ERROR_LOG)

log_msg(length(SNV_list), " cluster(s) with somatic SNVs")

# ── Répartition des clusters par étape de sortie ─────────────────────────────
if (file.exists(STAGE_LOG) && file.size(STAGE_LOG) > 0) {
    stage_df  <- fread(STAGE_LOG, header = FALSE, col.names = c("cluster_ID", "stage"))
    stage_tbl <- stage_df[, .N, by = stage][order(-N)]
    log_msg("Repartition des clusters par sortie :")
    for (i in seq_len(nrow(stage_tbl)))
        log_msg("  ", stage_tbl$stage[i], ": ", stage_tbl$N[i])
}


# =============================================================================
# 5. Map normal contigs to reference & report genomic coordinates
# =============================================================================

log_msg("Mapping contigs to reference genome...")

contig_fa  <- file.path(WD, "rawdata/contigs/contigs_normal.fa")
contig_sam <- file.path(WD, "rawdata/contigs/contigs_normal_aln.sam")
dir.create(dirname(contig_fa), showWarnings = FALSE, recursive = TRUE)

# Always write contig_normal[1] (forward strand) to the FASTA
writeLines(
    unlist(lapply(SNV_list, function(x)
        c(paste0(">contig_cluster_", x$cluster_ID), x$contig_normal[1])
    )),
    contig_fa
)

system(paste(MINIMAP, "-ax sr", REF, contig_fa, ">", contig_sam))

sam <- read.table(contig_sam, sep = "\t", comment.char = "@", fill = TRUE, row.names = NULL, stringsAsFactors = FALSE)

# ── Build SNV table ──────────────────────────────────────────────────────────

snv_rows <- lapply(seq_len(nrow(sam)), function(i) {
    flag <- as.integer(sam[i, 2])

    # Skip unmapped (0x4), secondary (0x100) and supplementary (0x800)
    # alignments. Secondary alignments repeat the same contig/mismatches at
    # a different locus and would otherwise duplicate SNV rows.
    if (is.na(flag) || bitwAnd(flag, 4) != 0 ||
        bitwAnd(flag, 256) != 0 || bitwAnd(flag, 2048) != 0)
        return(NULL)

    qname      <- as.character(sam[i, 1])
    cluster_ID <- as.integer(sub("contig_cluster_", "", qname))
    idx <- which(vapply(SNV_list, function(s) s$cluster_ID == cluster_ID, logical(1)))
    if (length(idx) == 0) return(NULL)
    x  <- SNV_list[[idx[1]]]
    mm <- x$mismatches

    is_reverse <- bitwAnd(flag, 16) != 0
    chrom      <- as.character(sam[i, 3])
    mapq       <- as.integer(sam[i, 5])
    ref_start  <- as.integer(sam[i, 4])
    cigar      <- as.character(sam[i, 6])
    contig_len <- nchar(x$contig_normal[1])

    need_rc <- xor(x$normal_rc, is_reverse)

    lapply(seq_len(nrow(mm)), function(j) {

        pat_pos <- mm$PatternStart[j]
        contig_fwd_pos <- if (x$normal_rc) contig_len - pat_pos + 1L else pat_pos
        query_pos <- if (is_reverse) contig_len - contig_fwd_pos + 1L else contig_fwd_pos

        ref_pos <- query_pos_to_ref_pos(query_pos, cigar, ref_start)
        if (is.na(ref_pos)) return(NULL)

        ref_nt <- as.character(mm$PatternSubstring[j])
        alt_nt <- as.character(mm$SubjectSubstring[j])
        if (need_rc) {
            ref_nt <- as.character(reverseComplement(DNAString(ref_nt)))
            alt_nt <- as.character(reverseComplement(DNAString(alt_nt)))
        }

        vaf_row <- x$vaf_df[j, ]

        data.frame(
            CHROM     = chrom,
            POS       = ref_pos,
            REF       = ref_nt,
            ALT       = alt_nt,
            n_alt            = vaf_row$n_alt,
            n_ref            = vaf_row$n_ref,
            VAF              = round(vaf_row$vaf, 4),
            vaf_p_value      = vaf_row$p_value,
            n_sig_kmers_in_normal = x$n_sig_kmers_in_normal,
            aln_score = x$aln_score,
            aln_score_norm = x$aln_score_norm,
            pid_aln   = x$pid_aln,
            n_kmer_sig = x$n_kmer_sig,
            n_kmer_sig_dist_k = x$n_kmer_sig_dist_k,
            n_reads   = x$n_reads,
            n_normal_reads = x$n_normal_reads,
            n_normal_kmers_total  = x$n_normal_kmers_total,
            n_normal_kmers_cov_ok = x$n_normal_kmers_cov_ok,
            frac_normal_kmers_cov_ok = x$frac_normal_kmers_cov_ok,
            n_mismatches_total   = x$n_mismatches_total,
            n_mismatches_somatic = x$n_mismatches_somatic,
            mean_weight = x$mean_weight,
            length_contig_tumour = x$length_contig_tumour,
            length_contig_normal = x$length_contig_normal,
            mapq        = mapq,
            cluster   = cluster_ID,
            stringsAsFactors = FALSE
        )
    })
})

snv_flat <- Filter(Negate(is.null), unlist(snv_rows, recursive = FALSE))

if (length(snv_flat) > 0) {
    snv_df <- do.call(rbind, snv_flat)
    snv_df <- snv_df[order(snv_df$CHROM, snv_df$POS), ]
} else {
    snv_df <- data.frame(
        CHROM = character(), POS = integer(),
        REF = character(), ALT = character(),
        n_alt = integer(), n_ref = integer(), VAF = numeric(),
        vaf_p_value = numeric(),
        n_sig_kmers_in_normal = numeric(),
        aln_score = numeric(), aln_score_norm = numeric(), pid_aln = numeric(),
        n_kmer_sig = numeric(), n_kmer_sig_dist_k = numeric(),
        n_reads = numeric(), n_normal_reads = numeric(),
        n_normal_kmers_total = numeric(), n_normal_kmers_cov_ok = numeric(),
        frac_normal_kmers_cov_ok = numeric(),
        n_mismatches_total = numeric(), n_mismatches_somatic = numeric(),
        mean_weight = numeric(),
        length_contig_tumour = numeric(), length_contig_normal = numeric(),
        mapq = integer(),
        cluster = integer(),
        stringsAsFactors = FALSE
    )
    log_msg("WARNING: no SNVs survived coordinate mapping")
}

# ── Write filtered TSV (final) ───────────────────────────────────────────────
header <- paste(c("#CHROM", "POS", "REF", "ALT",
                  "n_alt", "n_ref", "VAF", "vaf_p_value",
                  "n_sig_kmers_in_normal",
                  "aln_score", "aln_score_norm", "pid_aln",
                  "n_kmer_tumour_sig", "n_kmer_sig_dist_k",
                  "n_tumour_reads", "n_normal_reads",
                  "n_normal_kmers_total", "n_normal_kmers_cov_ok",
                  "frac_normal_kmers_cov_ok",
                  "n_mismatches_total", "n_mismatches_somatic",
                  "mean_weight",
                  "length_contig_tumour", "length_contig_normal",
                  "mapq",
                  "cluster_ID"), collapse = "\t")
tsv_out <- file.path(WD,
    paste0("output/results/somatic_snvs_DBG_", label_file, "_stats.tsv"))
writeLines(header, tsv_out)
write.table(snv_df, tsv_out, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE, append = TRUE)
log_msg(nrow(snv_df), " SNV(s) written to ", tsv_out, " (filtered)")