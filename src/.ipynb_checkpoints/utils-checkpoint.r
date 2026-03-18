# =============================================================================
# utils.R
# Helper functions for k-mer extraction, contig assembly, and alignment.
# Sourced by pipeline.R.
# =============================================================================

library(Biostrings)
library(pwalign)

WD        <- "/srv/home/mlef0011/VDARK"
KMC       <- file.path(WD, "software/kmc/bin/kmc")
KMC_TOOLS <- file.path(WD, "software/kmc/bin/kmc_tools")
TMP       <- file.path(WD, "tmp")

log_msg <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)


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