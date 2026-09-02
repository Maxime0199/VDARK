# =============================================================================
# assemble_dbg_igraph.r  —  De Bruijn Graph assembler (igraph engine)
#
# Drop-in replacement for assemble_dbg.r. Same public interface:
#   assemble_dbg(kmers, k, kmer_freq, coverage, mode) ->
#       list(contigs, bubbles, stats)
#
# Why this rewrite: the original environment/hashmap engine had a class of
# NULL vs character(0) bugs in .remove_tips() (setdiff(NULL, x) returns NULL
# in base R, which is indistinguishable from "found a real branch" -- see
# post-mortem). igraph's degree()/incident()/neighbors() always return
# length-0 vectors, never NULL, so that whole bug class is structurally
# impossible here. The graph-walk LOGIC (tip criteria, bubble criteria,
# unitig-extension criteria) is unchanged from assemble_dbg.r -- only the
# underlying graph representation changed.
# =============================================================================

suppressPackageStartupMessages(library(igraph))
# Biostrings/pwalign are already loaded by utils.r, sourced before this file.

# ── Helpers ───────────────────────────────────────────────────────────────────
.rc_one    <- function(s) as.character(reverseComplement(DNAString(s)))
.canonical <- function(km) { rc <- .rc_one(km); if (km < rc) km else rc }

# ── 1. Graph construction ─────────────────────────────────────────────────────
# Vertices = (k-1)-mers. Edges = k-mers (both strands), directed prefix->suffix.
# Edge attributes: kmer (sequence), freq (coverage).
.build_dbg_igraph <- function(kmers, k, kmer_freq) {

    canon_freq <- setNames(integer(0), character(0))
    if (!is.null(kmer_freq) && length(kmer_freq) > 0L) {
        nm <- vapply(names(kmer_freq), .canonical, character(1))
        canon_freq <- setNames(as.integer(kmer_freq), nm)
    }

    rc_kmers <- as.character(reverseComplement(DNAStringSet(kmers)))
    all_km   <- unique(c(kmers, rc_kmers))          # both strands, deduplicated

    pfx  <- substring(all_km, 1L, k - 1L)
    sfx  <- substring(all_km, 2L, k)
    freq <- vapply(all_km, function(km) {
        cf <- canon_freq[.canonical(km)]
        if (is.na(cf)) 1L else as.integer(cf)
    }, integer(1))

    edges_df <- data.frame(from = pfx, to = sfx, kmer = all_km, freq = freq,
                           stringsAsFactors = FALSE)

    g <- graph_from_data_frame(edges_df, directed = TRUE)
    E(g)$kmer <- edges_df$kmer
    E(g)$freq <- edges_df$freq
    g
}

# ── 2. Tip removal ────────────────────────────────────────────────────────────
# A tip = short (<= max_len edges), low-coverage (mean freq < cov_threshold)
# simple path that dead-ends (degree-1 terminus) and connects to the rest of
# the graph via a genuine branch (degree >= 2). Reaching the true, isolated
# end of a component (no branch anywhere) is NOT a tip and is left alone --
# this is exactly the distinction the original engine got wrong.
.trace_tip_from_sink <- function(g, sink_v, max_len, cov_threshold) {
    in_e <- as.integer(incident(g, sink_v, mode = "in"))
    for (e0 in in_e) {
        path  <- e0
        freqs <- E(g)$freq[e0]
        cur   <- as.integer(tail_of(g, e0))
        valid <- FALSE
        while (length(path) <= max_len) {
            if (degree(g, cur, mode = "out") >= 2) { valid <- TRUE; break }
            preds <- as.integer(incident(g, cur, mode = "in"))
            if (length(preds) != 1L) break     # 0 = true start (not a tip); >1 = merge (leave alone)
            e_next <- preds[1]
            path   <- c(e_next, path)
            freqs  <- c(E(g)$freq[e_next], freqs)
            cur    <- as.integer(tail_of(g, e_next))
        }
        if (valid && length(path) <= max_len && mean(freqs) < cov_threshold)
            return(path)
    }
    NULL
}

.trace_tip_from_source <- function(g, src_v, max_len, cov_threshold) {
    out_e <- as.integer(incident(g, src_v, mode = "out"))
    for (e0 in out_e) {
        path  <- e0
        freqs <- E(g)$freq[e0]
        cur   <- as.integer(head_of(g, e0))
        valid <- FALSE
        while (length(path) <= max_len) {
            if (degree(g, cur, mode = "in") >= 2) { valid <- TRUE; break }
            succs <- as.integer(incident(g, cur, mode = "out"))
            if (length(succs) != 1L) break
            e_next <- succs[1]
            path   <- c(path, e_next)
            freqs  <- c(freqs, E(g)$freq[e_next])
            cur    <- as.integer(head_of(g, e_next))
        }
        if (valid && length(path) <= max_len && mean(freqs) < cov_threshold)
            return(path)
    }
    NULL
}

#' Remove tips one at a time, rescanning after every removal.
#' Rescanning (rather than batching) sidesteps edge-ID renumbering after
#' delete_edges() -- small per-cluster graphs make this cheap.
.remove_tips_igraph <- function(g, max_len, cov_threshold) {
    n_removed <- 0L
    repeat {
        outdeg <- degree(g, mode = "out")
        indeg  <- degree(g, mode = "in")
        sinks   <- V(g)[outdeg == 0 & indeg > 0]
        sources <- V(g)[indeg == 0 & outdeg > 0]

        removed_path <- NULL
        for (v in sinks) {
            removed_path <- .trace_tip_from_sink(g, v, max_len, cov_threshold)
            if (!is.null(removed_path)) break
        }
        if (is.null(removed_path)) {
            for (v in sources) {
                removed_path <- .trace_tip_from_source(g, v, max_len, cov_threshold)
                if (!is.null(removed_path)) break
            }
        }
        if (is.null(removed_path)) break

        g <- delete_edges(g, removed_path)
        n_removed <- n_removed + length(removed_path)
    }
    list(graph = g, n_removed = n_removed)
}

# ── 3. Bubble popping ─────────────────────────────────────────────────────────
# Bubble = two paths from a branch node (out-degree >= 2) converging at a
# common downstream node within max_depth. Majority-coverage path kept;
# minority recorded (for germline pre-filtering) and deleted.
.path_to_seq_igraph <- function(g, edge_ids, k) {
    if (length(edge_ids) == 0L) return("")
    kmers <- E(g)$kmer[edge_ids]
    paste0(kmers[1L], paste(substring(kmers[-1L], k, k), collapse = ""))
}

.pop_bubbles_igraph <- function(g, max_depth, k) {
    records <- list()
    branch_nodes <- as.integer(V(g)[degree(g, mode = "out") >= 2])

    for (B in branch_nodes) {
        out_e <- as.integer(incident(g, B, mode = "out"))
        if (length(out_e) < 2L) next   # may already have changed this call

        explore <- function(e0) {
            edges <- e0
            cur   <- as.integer(head_of(g, e0))
            nodes <- cur
            for (d in seq_len(max_depth - 1L)) {
                succ <- as.integer(incident(g, cur, mode = "out"))
                if (length(succ) != 1L) break
                e_next <- succ[1]
                edges  <- c(edges, e_next)
                cur    <- as.integer(head_of(g, e_next))
                nodes  <- c(nodes, cur)
            }
            list(edges = edges, nodes = nodes)
        }

        paths      <- lapply(out_e, explore)
        node_lists <- lapply(paths, `[[`, "nodes")
        all_nodes  <- unique(unlist(node_lists))
        common     <- all_nodes[vapply(all_nodes, function(n)
            sum(vapply(node_lists, function(nl) n %in% nl, logical(1))) >= 2L,
            logical(1))]
        if (length(common) == 0L) next

        depths <- vapply(common, function(n)
            max(vapply(node_lists, function(nl) {
                i <- match(n, nl); if (is.na(i)) Inf else i
            }, numeric(1))), numeric(1))
        merge_node <- common[which.min(depths)]

        truncated <- Filter(Negate(is.null), lapply(paths, function(p) {
            j <- match(merge_node, p$nodes)
            if (is.na(j)) return(NULL)
            e <- p$edges[seq_len(j)]
            list(edges = e, freqs = E(g)$freq[e])
        }))
        if (length(truncated) < 2L) next

        totals   <- vapply(truncated, function(p) sum(p$freqs), numeric(1))
        keep_idx <- which.max(totals)
        kept     <- truncated[[keep_idx]]

        for (i in seq_along(truncated)) {
            if (i == keep_idx) next
            rem <- truncated[[i]]
            if (length(rem$edges) != length(kept$edges)) next  # indel bubble: skip (SNP-only, as original)
            records[[length(records) + 1L]] <- list(
                kept_seq     = .path_to_seq_igraph(g, kept$edges, k),
                removed_seq  = .path_to_seq_igraph(g, rem$edges,  k),
                kept_freq    = mean(kept$freqs),
                removed_freq = mean(rem$freqs)
            )
            g <- delete_edges(g, rem$edges)
        }
    }
    list(graph = g, n_bubbles = length(records), records = records)
}

# ── 4. Unitig extraction ──────────────────────────────────────────────────────
# Walk linearly from a seed edge; stop at any unresolved branch or merge point.
.extract_unitig_igraph <- function(g, seed_edge, k) {
    visited_e <- seed_edge

    fwd <- seed_edge
    cur <- as.integer(head_of(g, seed_edge))
    repeat {
        succ <- as.integer(incident(g, cur, mode = "out"))
        succ <- succ[!(succ %in% visited_e)]
        if (length(succ) != 1L) break
        nxt <- succ[1]
        nxt_target <- as.integer(head_of(g, nxt))
        if (degree(g, nxt_target, mode = "in") > 1L) break
        fwd <- c(fwd, nxt)
        visited_e <- c(visited_e, nxt)
        cur <- nxt_target
    }

    bwd <- integer(0)
    cur <- as.integer(tail_of(g, seed_edge))
    repeat {
        pred <- as.integer(incident(g, cur, mode = "in"))
        pred <- pred[!(pred %in% visited_e)]
        if (length(pred) != 1L) break
        prv <- pred[1]
        prv_source <- as.integer(tail_of(g, prv))
        if (degree(g, prv_source, mode = "out") > 1L) break
        bwd <- c(prv, bwd)
        visited_e <- c(visited_e, prv)
        cur <- prv_source
    }

    .path_to_seq_igraph(g, c(bwd, fwd), k)
}

# ── 5. Bubble alleles → contig coordinates (unchanged from assemble_dbg.r) ───
.empty_bubbles_df <- function()
    data.frame(pos_in_contig = integer(), ref_allele = character(),
               alt_allele = character(), ref_freq = numeric(),
               alt_freq = numeric(), stringsAsFactors = FALSE)

.bubbles_to_snps <- function(records, contig_fwd) {
    if (length(records) == 0L) return(.empty_bubbles_df())
    contig_rc <- .rc_one(contig_fwd)
    rows      <- list()
    for (b in records) {
        kept_seq    <- b$kept_seq
        removed_seq <- b$removed_seq
        if (nchar(kept_seq) != nchar(removed_seq)) next
        diffs <- which(strsplit(kept_seq,    "")[[1L]] !=
                       strsplit(removed_seq, "")[[1L]])
        if (length(diffs) == 0L) next
        m_fwd <- regexpr(kept_seq, contig_fwd, fixed = TRUE)
        m_rc  <- regexpr(kept_seq, contig_rc,  fixed = TRUE)
        for (d in diffs) {
            ref_n <- substr(kept_seq,    d, d)
            alt_n <- substr(removed_seq, d, d)
            if (m_fwd[1L] > 0L) {
                rows[[length(rows) + 1L]] <- data.frame(
                    pos_in_contig = m_fwd[1L] + d - 1L,
                    ref_allele = ref_n, alt_allele = alt_n,
                    ref_freq = b$kept_freq, alt_freq = b$removed_freq,
                    stringsAsFactors = FALSE)
            } else if (m_rc[1L] > 0L) {
                rows[[length(rows) + 1L]] <- data.frame(
                    pos_in_contig = nchar(contig_fwd) - (m_rc[1L] + d - 1L) + 1L,
                    ref_allele = .rc_one(ref_n), alt_allele = .rc_one(alt_n),
                    ref_freq = b$kept_freq, alt_freq = b$removed_freq,
                    stringsAsFactors = FALSE)
            }
        }
    }
    if (length(rows) == 0L) return(.empty_bubbles_df())
    do.call(rbind, rows)
}

# ── Entry point (same signature as assemble_dbg.r) ────────────────────────────
assemble_dbg <- function(kmers, k = K, kmer_freq,
                         coverage, mode = c("normal", "tumour")) {
    mode <- match.arg(mode)
    if (length(kmers) == 0L) return(NULL)
    if (length(kmers) == 1L) return(list(
        contigs = c(kmers, .rc_one(kmers)),
        bubbles = .empty_bubbles_df(),
        stats   = list(n_kmers_input = 1L, n_tips_removed = 0L,
                       n_bubbles = 0L, contig_len = nchar(kmers))
    ))

    tip_cov      <- max(2, coverage / if (mode == "normal") 4 else 8)
    bubble_depth <- 2L * k

    g <- .build_dbg_igraph(kmers, k, kmer_freq)

    tip_res <- .remove_tips_igraph(g, max_len = 2L * k, cov_threshold = tip_cov)
    g <- tip_res$graph

    bub_res <- .pop_bubbles_igraph(g, max_depth = bubble_depth, k = k)
    g <- bub_res$graph

    if (ecount(g) == 0L) return(NULL)

    freqs      <- E(g)$freq
    seed_order <- order(freqs, decreasing = TRUE)[seq_len(min(5L, ecount(g)))]

    contig_fwd <- NA_character_
    for (s in seed_order) {
        cand <- .extract_unitig_igraph(g, s, k)
        if (!is.na(cand) && (is.na(contig_fwd) || nchar(cand) > nchar(contig_fwd)))
            contig_fwd <- cand
    }
    if (is.na(contig_fwd) || nchar(contig_fwd) < k) return(NULL)

    list(
        contigs = c(contig_fwd, .rc_one(contig_fwd)),
        bubbles = .bubbles_to_snps(bub_res$records, contig_fwd),
        stats   = list(n_kmers_input  = length(kmers),
                       n_tips_removed = tip_res$n_removed,
                       n_bubbles      = bub_res$n_bubbles,
                       contig_len     = nchar(contig_fwd))
    )
}