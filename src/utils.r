# =============================================================================
# utils.R  (corrected)
# Helper functions for k-mer extraction, contig assembly, alignment,
# CIGAR parsing, and VAF estimation.
# utils.R
# Helper functions for k-mer extraction, contig assembly, and alignment.
# Sourced by pipeline.R.
# =============================================================================

library(Biostrings)
library(pwalign)

K <- 31
WD        <- "/srv/home/mlef0011/VDARK"
KMC       <- file.path(WD, "software/kmc/bin/kmc")
KMC_TOOLS <- file.path(WD, "software/kmc/bin/kmc_tools")
TMP       <- file.path(WD, "tmp")
SPADES <- "/srv/home/mlef0011/anaconda3/envs/VDARK/bin/spades.py"


log_msg <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)


inspect_and_split_cluster <- function(cluster_ID, G, clusters,
                                      min_modularity = 0.25, resolution = 1.0) {

    reads_in <- names(clusters$membership)[clusters$membership == cluster_ID]
    subG     <- induced_subgraph(G, vids = reads_in)

    ap <- articulation_points(subG)
    cat("Reads:", vcount(subG), "| Edges:", ecount(subG),
        "| Points d'articulation:", length(ap), "\n")
    if (length(ap) > 0) print(V(subG)$name[ap])

    comm <- cluster_leiden(subG, objective_function = "modularity",
                           weights = E(subG)$weight, resolution = resolution)

    n_comm <- length(unique(membership(comm)))
    mod    <- if (n_comm > 1) modularity(subG, membership(comm), weights = E(subG)$weight) else NA_real_

    cat("Communautés Leiden:", n_comm, "| modularité:", round(mod, 3), "\n")

    if (n_comm <= 1 || is.na(mod) || mod < min_modularity) {
        cat("-> Pas de split fiable, cluster conservé tel quel.\n")
        return(invisible(list(split = FALSE)))
    }

    sub_clusters <- split(V(subG)$name, membership(comm))
    cat("-> Split en", length(sub_clusters), "sous-cluster(s):\n")
    for (i in seq_along(sub_clusters))
        cat("  sous-cluster", i, ":", length(sub_clusters[[i]]), "reads\n")

    invisible(list(split = TRUE, sub_clusters = sub_clusters,
                   modularity = mod, articulation_points = V(subG)$name[ap]))
}


# ── K-mer utilities ───────────────────────────────────────────────────────────

#' Extract canonical k-mers from a sequence
get_kmers <- function(seq, k) {
    n <- nchar(seq) - k + 1
    if (n <= 0) return(character(0))
    kmers <- substring(seq, 1:n, k:nchar(seq))
    rc    <- as.character(reverseComplement(DNAStringSet(kmers)))
    ifelse(kmers < rc, kmers, rc)
}


# ── Contig assembly ───────────────────────────────────────────────────────────

#' Greedy frequency-weighted k-mer assembly
assemble_kmers <- function(kmers, k, kmer_freq = NULL) {
    kmers <- as.character(kmers) 
    if (length(kmers) == 0) return(NULL)
    if (length(kmers) == 1)
        return(c(kmers, as.character(reverseComplement(DNAString(kmers)))))

    # Build frequency lookup
    freq_env <- new.env(hash = TRUE)
    if (is.null(kmer_freq)) {
        for (km in kmers) freq_env[[km]] <- 1L
    } else {
        for (i in seq_along(kmer_freq)) freq_env[[names(kmer_freq)[i]]] <- kmer_freq[i]
    }

    get_freq <- function(km) {
        f <- freq_env[[km]]
        if (!is.null(f)) return(f)
        rc <- as.character(reverseComplement(DNAString(km)))
        f  <- freq_env[[rc]]
        if (!is.null(f)) return(f)
        1L
    }

    # Expand to both strands
    rc_kmers  <- as.character(reverseComplement(DNAStringSet(kmers)))
    all_forms <- unique(c(kmers, rc_kmers))

    # Build suffix -> next k-mer(s) lookup
    suffix_map <- new.env(hash = TRUE)
    for (km in all_forms) {
        prefix   <- substring(km, 1, k - 1)
        existing <- suffix_map[[prefix]]
        suffix_map[[prefix]] <- if (is.null(existing)) km else c(existing, km)
    }

    # Identify start k-mers (prefix not found as suffix)
    all_prefixes <- substring(all_forms, 1, k - 1)
    all_suffixes <- substring(all_forms, 2, k)
    start_kmers  <- all_forms[!(all_prefixes %in% all_suffixes)]
    if (length(start_kmers) == 0) {
        start_kmers <- all_forms[1]
        log_msg("Cycle detected in k-mer graph")
    }

    # Greedy extension from each start k-mer
    contigs <- lapply(start_kmers, function(start) {
        contig  <- start
        current <- start
        visited <- new.env(hash = TRUE)
        visited[[current]] <- TRUE

        repeat {
            candidates <- suffix_map[[substring(current, 2, k)]]
            if (is.null(candidates)) break
            candidates <- candidates[vapply(candidates,
                function(km) is.null(visited[[km]]), logical(1))]
            if (length(candidates) == 0) break

            next_kmer <- if (length(candidates) == 1) candidates else {
                freqs <- vapply(candidates, get_freq, numeric(1))
                candidates[which.max(freqs)]
            }
            contig  <- paste0(contig, substring(next_kmer, k, k))
            visited[[next_kmer]] <- TRUE
            current <- next_kmer
        }
        contig
    })

    best    <- contigs[[which.max(nchar(unlist(contigs)))]]
    rc_best <- as.character(reverseComplement(DNAString(best)))
    c(best, rc_best)
}

                                     
                                     
# ── K-mer utilities ───────────────────────────────────────────────────────────

#' Extract canonical k-mers from a sequence
get_kmers <- function(seq, k = K) {
    n <- nchar(seq) - k + 1
    if (n <= 0) return(character(0))
    kmers <- substring(seq, 1:n, k:nchar(seq))
    rc    <- as.character(reverseComplement(DNAStringSet(kmers)))
    ifelse(kmers < rc, kmers, rc)
}

#' Extract ALL k-mers (non-canonical, both strands not collapsed)
#' Useful for VAF counting where we need exact matching
get_kmers_raw <- function(seq, k = K) {
    n <- nchar(seq) - k + 1
    if (n <= 0) return(character(0))
    substring(seq, 1:n, k:nchar(seq))
}


# ── CIGAR parsing & coordinate conversion ─────────────────────────────────────

#' Parse a CIGAR string into a data.frame of (length, operation) pairs
#' @return data.frame with columns: len (integer), op (character)
parse_cigar <- function(cigar) {
    ops <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))[[1]]
    data.frame(
        len = as.integer(sub("[A-Z=]+$", "", ops)),
        op  = sub("^[0-9]+", "", ops),
        stringsAsFactors = FALSE
    )
}

#' Convert a 1-based query position to a 1-based reference position using CIGAR.
#'
#' Walks through the CIGAR operations and tracks how query positions map
#' to reference positions.  Handles M/=/X (consume both), I (query only),
#' D/N (reference only), S (soft-clip, query only), H (hard-clip, neither).
#'
#' @param query_pos  1-based position in the query (= contig as in SAM SEQ)
#' @param cigar      CIGAR string from the SAM record
#' @param ref_start  POS field from the SAM record (1-based leftmost ref position)
#' @return integer reference position, or NA if query_pos falls in an
#'         insertion or clipped region
query_pos_to_ref_pos <- function(query_pos, cigar, ref_start) {
    parsed <- parse_cigar(cigar)

    q <- 0L              # query bases consumed so far (0-based counter)
    r <- ref_start - 1L  # reference position tracker  (0-based)

    for (i in seq_len(nrow(parsed))) {
        op  <- parsed$op[i]
        len <- parsed$len[i]

        if (op %in% c("M", "=", "X")) {
            # Consumes both query and reference
            if (query_pos <= q + len) {
                return(r + (query_pos - q))
            }
            q <- q + len
            r <- r + len

        } else if (op == "I") {
            # Insertion — consumes query only
            if (query_pos <= q + len) return(NA_integer_)
            q <- q + len

        } else if (op == "S") {
            # Soft clip — consumes query only (bases present in SEQ but unaligned)
            if (query_pos <= q + len) return(NA_integer_)
            q <- q + len

        } else if (op %in% c("D", "N")) {
            # Deletion / skipped region — consumes reference only
            r <- r + len

        }
        # H (hard clip) and P (padding) consume neither → skip
    }
    NA_integer_
}



                                            

# ── Germline filter ───────────────────────────────────────────────────────────


is_germline <- function(mm, contig_strand, normal_seqs, k = K, threshold = 10,
                        buf = log_msg) {

    if (nrow(mm) == 0) return(logical(0))

    reads_dna <- DNAStringSet(normal_seqs)

    # ── Groupes de mismatches dont les fenêtres de k-mers se chevauchent ────
    # (positions distantes de moins de k pb -> même fenêtre potentielle)
    ord      <- order(mm$PatternStart)
    pos_ord  <- mm$PatternStart[ord]
    new_grp  <- c(TRUE, diff(pos_ord) >= k)
    group_id <- integer(nrow(mm))
    group_id[ord] <- cumsum(new_grp)

    count_matching_reads <- function(seq, pos, k) {
        s_min <- max(1, pos - k + 1)
        s_max <- min(pos, nchar(seq) - k + 1)
        if (s_min > s_max) return(0L)
        mut_kmers <- substring(seq, s_min:s_max, (s_min:s_max) + k - 1)
        rc_mut    <- as.character(reverseComplement(DNAStringSet(mut_kmers)))
        pdict     <- PDict(unique(c(mut_kmers, rc_mut)))
        sum(colSums(vcountPDict(pdict, reads_dna)) > threshold)
    }

    vapply(seq_len(nrow(mm)), function(i) {

        pos <- mm$PatternStart[i]
        if (pos < 1 || pos > nchar(contig_strand)) return(FALSE)
        if (substr(contig_strand, pos, pos) != as.character(mm$PatternSubstring[i]))
            return(FALSE)

        # ── Test A : seule cette position mutée, voisins laissés tels quels ──
        seq_a <- contig_strand
        substr(seq_a, pos, pos) <- as.character(mm$SubjectSubstring[i])
        n_reads_a <- count_matching_reads(seq_a, pos, k)

        # ── Test B : cette position + tous ses voisins de groupe mutés ensemble ──
        mates <- which(group_id == group_id[i])
        n_reads_b <- 0L
        if (length(mates) > 1) {
            seq_b <- contig_strand
            for (m in mates)
                substr(seq_b, mm$PatternStart[m], mm$PatternStart[m]) <-
                    as.character(mm$SubjectSubstring[m])
            n_reads_b <- count_matching_reads(seq_b, pos, k)
        }

        n_reads <- max(n_reads_a, n_reads_b)
        is_germ <- n_reads > 0
        buf("  pos ", pos, " | A=", n_reads_a, " B=", n_reads_b,
            " normal reads | ", if (is_germ) "germline" else "somatic")
        is_germ
    }, logical(1))
}
                                            
# ── Alignment ─────────────────────────────────────────────────────────────────

#' Best local alignment across all 4 strand orientations.
#'
#' Returns a list with:
#'   $alignment   — the PairwiseAlignmentsSingleSubject object
#'   $normal_rc   — TRUE if contig_normal[2] (RC) was the best pattern
#'
#' Scoring rationale:
#'   match=1, mismatch=-1  — standard for closely related sequences (tumour
#'     vs normal differ by few SNVs), balanced so that a single SNV doesn't
#'     mask surrounding matches.
#'   gapOpening=-5, gapExtension=-2 — penalises gaps (indels in the contig
#'     alignment) enough to avoid spurious gapped alignments for SNV-only
#'     detection.  If you ever extend to indel calling, consider relaxing
#'     gapOpening to -3.
best_alignment <- function(contig_normal, contig_tumour) {
    orientations <- list(
        list(n = contig_normal[1], t = contig_tumour[1], normal_rc = FALSE),
        list(n = contig_normal[1], t = contig_tumour[2], normal_rc = FALSE),
        list(n = contig_normal[2], t = contig_tumour[1], normal_rc = TRUE),
        list(n = contig_normal[2], t = contig_tumour[2], normal_rc = TRUE)
    )
    alns <- lapply(orientations, function(o) {
        pairwiseAlignment(
            DNAString(o$n), DNAString(o$t),


# ── Germline filter ───────────────────────────────────────────────────────────

#' Returns TRUE for each mismatch found in the normal sample (germline).
#' buf: logging function from the caller (e.g. buf <- function(...) logs <<- c(logs, ...))
is_germline <- function(mm, contig_normal, NORMAL_R1, NORMAL_R2, cluster_ID, k , threshold = 3,
                        buf = log_msg) {

    dir.create(file.path(TMP, cluster_ID), showWarnings = FALSE)

    vapply(seq_len(nrow(mm)), function(i) {
        pos <- mm$PatternStart[i]
        seq <- contig_normal[1]

        if (substr(seq, pos, pos) != as.character(mm$PatternSubstring[i]))
            return(FALSE)

        substr(seq, pos, pos) <- as.character(mm$SubjectSubstring[i])

        s_min     <- max(1, pos - k + 1)
        s_max     <- min(pos, nchar(seq) - k + 1)
        mut_kmers <- substring(seq, s_min:s_max, (s_min:s_max) + k - 1)

        tag      <- paste0("cl", cluster_ID, "_pos", pos)
        fa_file  <- file.path(TMP, cluster_ID, paste0("mut_kmers_", tag, ".fa"))
        kmc_db_g <- file.path(TMP, cluster_ID, paste0("mut_kmc_", tag))
        fq_R1_g  <- file.path(TMP, cluster_ID, paste0("normal_R1_", tag, ".fq"))
        fq_R2_g  <- file.path(TMP, cluster_ID, paste0("normal_R2_", tag, ".fq"))

        writeLines(paste0(">k_", seq_along(mut_kmers), "\n", mut_kmers), fa_file)
        system(paste(KMC, paste0("-k",K),"-t12 -ci1 -fm", fa_file, kmc_db_g, TMP),
               ignore.stdout = TRUE, ignore.stderr = TRUE)
        system(paste(KMC_TOOLS, "filter", kmc_db_g, "-ci1", NORMAL_R1, "-ci1", fq_R1_g),
               ignore.stdout = TRUE, ignore.stderr = TRUE)
        system(paste(KMC_TOOLS, "filter", kmc_db_g, "-ci1", NORMAL_R2, "-ci1", fq_R2_g),
               ignore.stdout = TRUE, ignore.stderr = TRUE)

        reads <- tryCatch(
            as.character(c(readDNAStringSet(fq_R1_g, format = "fastq"),
                           readDNAStringSet(fq_R2_g, format = "fastq"))),
            error = function(e) character(0)
        )
        file.remove(fa_file, fq_R1_g, fq_R2_g)

        if (length(reads) == 0) return(FALSE)

        rc_mut  <- as.character(reverseComplement(DNAStringSet(mut_kmers)))
        pdict   <- PDict(unique(c(mut_kmers, rc_mut)))
        n_reads <- sum(colSums(vcountPDict(pdict, DNAStringSet(reads))) > 0)
        label   <- if (n_reads > threshold) "germline" else "somatic"
        buf("  pos ", pos, " | ", n_reads, " normal reads | ", label)
        
        n_reads > threshold
    }, logical(1))
}


# ── Alignment ─────────────────────────────────────────────────────────────────

#' Best local alignment across all 4 strand orientations
best_alignment <- function(contig_normal, contig_tumour) {
    orientations <- list(
        c(contig_normal[1], contig_tumour[1]),
        c(contig_normal[1], contig_tumour[2]),
        c(contig_normal[2], contig_tumour[1]),
        c(contig_normal[2], contig_tumour[2])
    )
    alns <- lapply(orientations, function(x) {
        pairwiseAlignment(
            DNAString(x[1]), DNAString(x[2]),
            type = "local",
            substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = -1),
            gapOpening   = -5,
            gapExtension = -2
        )
    })
    best_idx <- which.max(sapply(alns, score))
    list(
        alignment = alns[[best_idx]],
        normal_rc = orientations[[best_idx]]$normal_rc
    )

    alns[[which.max(sapply(alns, score))]]
}


#' Format alignment into buffer (or print directly if buf = log_msg)
format_alignment <- function(aln, cluster_ID, buf = log_msg) {
    pat  <- strsplit(as.character(pattern(aln)), "")[[1]]
    sub_ <- strsplit(as.character(subject(aln)), "")[[1]]
    diff <- ifelse(pat == sub_, ".", "*")
    buf("  --- Alignment cluster ", cluster_ID, " ---")
    buf("  Normal: ", paste(pat,  collapse = ""))
    buf("          ", paste(diff, collapse = ""))
    buf("  Tumour: ", paste(sub_, collapse = ""))
}


# ── VAF estimation ────────────────────────────────────────────────────────────

#' Extract ALL tumour reads at a locus using flanking k-mers.

extract_tumour_reads_at_locus <- function(kmc_db, tmp_dir,
                                          min_hits = 10) {
 
    TUMOUR_R1 <- file.path(WD, "rawdata/reads/tumour_chr21_R1.fq")
    TUMOUR_R2 <- file.path(WD, "rawdata/reads/tumour_chr21_R2.fq")
 
    fq_R1 <- file.path(tmp_dir, paste0("tumour_vaf_R1_.fq"))
    fq_R2 <- file.path(tmp_dir, paste0("tumour_vaf_R2_.fq"))
 
    ci_flag <- paste0("-ci", min_hits)
    system(paste(KMC_TOOLS, "filter", kmc_db, "-ci1", TUMOUR_R1, ci_flag, fq_R1),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste(KMC_TOOLS, "filter", kmc_db, "-ci1", TUMOUR_R2, ci_flag, fq_R2),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
 
    reads <- tryCatch(
        as.character(c(readDNAStringSet(fq_R1, format = "fastq"),
                       readDNAStringSet(fq_R2, format = "fastq"))),
        error = function(e) character(0)
    )
    file.remove(fq_R1, fq_R2)
    reads
}

#' Estimate variant allele frequency (VAF) from ALL tumour reads at a locus.
                                            
estimate_vaf <- function(mm, contig_strand, all_tumour_reads, k = K) {
 
    na_row <- data.frame(n_alt = NA_integer_, n_ref = NA_integer_, vaf = NA_real_)
 
    if (length(all_tumour_reads) == 0)
        return(do.call(rbind, replicate(nrow(mm), na_row, simplify = FALSE)))
 
    reads_dna <- DNAStringSet(all_tumour_reads)
 
    rows <- lapply(seq_len(nrow(mm)), function(i) {
 
        pos        <- mm$PatternStart[i]
        ref_allele <- as.character(mm$PatternSubstring[i])
        alt_allele <- as.character(mm$SubjectSubstring[i])
        seq_len    <- nchar(contig_strand)
 
        # ── Build REF and ALT k-mers from the FULL contig ────────────────────
        s_min <- max(1L, pos - k + 1L)
        s_max <- min(pos, seq_len - k + 1L)
        if (s_min > s_max) return(na_row)
 
        ref_kmers <- substring(contig_strand, s_min:s_max, (s_min:s_max) + k - 1L)
 
        alt_contig <- contig_strand
        substr(alt_contig, pos, pos) <- alt_allele
        alt_kmers <- substring(alt_contig, s_min:s_max, (s_min:s_max) + k - 1L)
 
        # Sanity: remove any k-mers with non-ACGT characters
        ref_kmers <- ref_kmers[grepl("^[ACGT]+$", ref_kmers)]
        alt_kmers <- alt_kmers[grepl("^[ACGT]+$", alt_kmers)]
 
        if (length(ref_kmers) == 0 || length(alt_kmers) == 0) return(na_row)
 
        # ── Include reverse complements ──────────────────────────────────────
        ref_all <- unique(c(ref_kmers,
                            as.character(reverseComplement(DNAStringSet(ref_kmers)))))
        alt_all <- unique(c(alt_kmers,
                            as.character(reverseComplement(DNAStringSet(alt_kmers)))))
 
        # ── Count per read ───────────────────────────────────────────────────
        tryCatch({
            pdict_ref <- PDict(DNAStringSet(ref_all))
            pdict_alt <- PDict(DNAStringSet(alt_all))
 
            ref_counts <- colSums(vcountPDict(pdict_ref, reads_dna))
            alt_counts <- colSums(vcountPDict(pdict_alt, reads_dna))
 
            informative <- (ref_counts + alt_counts) > 0
            if (sum(informative) == 0) return(na_row)
 
            # Assign each read to ALT or REF by majority of hits
            # Ties → REF (conservative)
            n_alt <- sum(alt_counts[informative] > ref_counts[informative])
            n_ref <- sum(ref_counts[informative] >= alt_counts[informative])
 
            data.frame(n_alt = n_alt, n_ref = n_ref,
                       vaf   = n_alt / (n_alt + n_ref))
        }, error = function(e) na_row)
    })
 
    do.call(rbind, rows)
}
 
                   

# Filter mismatches whose (pos, ref, alt) matches a normal-contig bubble.
filter_germline_bubbles <- function(mm, contig_strand, normal_bubbles_df,
                                    normal_rc, buf = log_msg) {

    if (is.null(normal_bubbles_df) || nrow(normal_bubbles_df) == 0L)
        return(rep(FALSE, nrow(mm)))

    contig_len    <- nchar(contig_strand)
    bubbles_local <- normal_bubbles_df

    if (normal_rc) {
        # Flip positions to RC strand
        bubbles_local$pos_in_contig <-
            contig_len - bubbles_local$pos_in_contig + 1L
        # Flip allele orientation
        bubbles_local$ref_allele <- as.character(reverseComplement(
            DNAStringSet(bubbles_local$ref_allele)))
        bubbles_local$alt_allele <- as.character(reverseComplement(
            DNAStringSet(bubbles_local$alt_allele)))
    }

    # Tolerate ±1 nt to absorb any off-by-one in path locating
    is_germ <- vapply(seq_len(nrow(mm)), function(i) {
        pos <- mm$PatternStart[i]
        ref <- as.character(mm$PatternSubstring[i])
        alt <- as.character(mm$SubjectSubstring[i])
        any(bubbles_local$pos_in_contig %in% (pos + (-1L:1L)) &
            bubbles_local$ref_allele == ref &
            bubbles_local$alt_allele == alt)
    }, logical(1))

    n_filtered <- sum(is_germ)
    if (n_filtered > 0L)
        buf("  Pre-filtered ", n_filtered,
            " mismatches matching normal-contig bubbles (germline)")
    is_germ
}

