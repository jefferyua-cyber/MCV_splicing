library(Rsubread)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)

# ---- Parameters -------------------------------------------------------

REFERENCE_FASTA <- ""

FASTQ_DIR    <- ""

SAMPLE_NAMES <- c("")

ALIGN_INDEX  <- ""

BAM_SUFFIX   <- "aligned"

PAIRED_END   <- TRUE

SEQ_TYPE     <- "rna"

OUTPUT_DIR   <- ""

SPLIT_STRAND <- TRUE

# ---- Step functions ----------------------------------------------------

run_step0_buildindex <- function(align_index, reference_fasta) {
  message("\n", strrep("=", 60))
  message("STEP 0: Building genome index")
  message(strrep("=", 60))
  if (!file.exists(reference_fasta))
    stop("FASTA not found: ", reference_fasta,
         "\n  Update REFERENCE_FASTA to the correct path.")
  Rsubread::buildindex(basename = align_index, reference = reference_fasta)
  message("Done. Comment out run_step0_buildindex() for future runs.")
}

run_step1_align <- function(fastq_dir, output_dir, align_index, bam_suffix,
                            sample_names = NULL, paired_end = FALSE,
                            seq_type = "rna") {
  message("\n", strrep("=", 60))
  message("STEP 1: Aligning FASTQ -> BAM")
  message(strrep("=", 60))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  fastq_r1 <- sort(list.files(fastq_dir, pattern = "_R1_001\\.fastq\\.gz$",
                              full.names = TRUE))
  if (length(fastq_r1) == 0)
    stop("No *_R1_001.fastq.gz files found in: ", fastq_dir,
         "\n  Check FASTQ_DIR and that filenames end in _R1_001.fastq.gz")

  snames <- if (!is.null(sample_names)) {
    if (length(sample_names) != length(fastq_r1))
      stop("SAMPLE_NAMES has ", length(sample_names), " entries but ",
           length(fastq_r1), " FASTQ files found (alphabetical order):\n",
           paste(" ", basename(fastq_r1), collapse = "\n"))
    sample_names
  } else {
    gsub("_R1_001\\.fastq\\.gz$", "", basename(fastq_r1))
  }

  bam_paths <- character(length(fastq_r1))
  for (i in seq_along(fastq_r1)) {
    r1  <- fastq_r1[i]
    bam <- file.path(output_dir, paste0(snames[i], "_", bam_suffix, ".bam"))
    bam_paths[i] <- bam
    message(sprintf("\n  [%d/%d] %s -> %s",
                    i, length(fastq_r1), basename(r1), basename(bam)))
    if (paired_end) {
      r2 <- gsub("_R1_001\\.fastq\\.gz$", "_R2_001.fastq.gz", r1)
      if (!file.exists(r2)) stop("R2 file not found: ", r2)
      Rsubread::align(index = align_index, readfile1 = r1, readfile2 = r2,
                      output_file = bam, type = seq_type)
    } else {
      Rsubread::align(index = align_index, readfile1 = r1,
                      output_file = bam, type = seq_type)
    }
    message("  Done: ", basename(bam))
  }
  message("\nStep 1 complete. BAMs saved to: ", output_dir)
  return(invisible(bam_paths))
}

run_step2_sort_index <- function(bam_paths) {
  message("\n", strrep("=", 60))
  message("STEP 2: Sorting and indexing BAMs")
  message(strrep("=", 60))
  for (bam in bam_paths) {
    if (file.exists(paste0(bam, ".bai"))) {
      message("  Already indexed, skipping: ", basename(bam))
      next
    }
    message("  Sorting:  ", basename(bam))
    sorted_prefix <- sub("\\.bam$", "_sorted", bam)
    Rsamtools::sortBam(bam, sorted_prefix)
    file.rename(paste0(sorted_prefix, ".bam"), bam)
    message("  Indexing: ", basename(bam))
    Rsamtools::indexBam(bam)
    message("  Done.")
  }
  message("\nStep 2 complete.")
  return(invisible(bam_paths))
}

run_step3_strandedness <- function(bamfile, n_reads = 200000, paired_end = FALSE) {
  message("\n", strrep("=", 60))
  message("STEP 3: Detecting library strandedness")
  message(strrep("=", 60))
  message("  BAM: ", basename(bamfile))

  flags <- Rsamtools::scanBam(bamfile, param = Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE),
    what = "flag"))[[1]]$flag

  if (length(flags) == 0)
    stop("No mapped reads in: ", basename(bamfile),
         "\n  Check that Step 1 completed successfully.")

  flags      <- flags[seq_len(min(n_reads, length(flags)))]
  is_reverse <- bitwAnd(flags, 16L) > 0
  is_read1   <- bitwAnd(flags, 64L) > 0

  if (!paired_end) {
    pct_fwd <- round(100 * sum(!is_reverse) / length(flags), 1)
    pct_rev <- round(100 * sum( is_reverse) / length(flags), 1)
    message(sprintf("  Forward reads : %.1f%%  |  Reverse reads : %.1f%%",
                    pct_fwd, pct_rev))
    result <- if (pct_fwd >= 70) "forward stranded" else
      if (pct_rev >= 70) "reverse stranded" else "unstranded"
  } else {
    r1_flags <- flags[is_read1]
    pct_fwd  <- round(100 * sum(bitwAnd(r1_flags, 16L) == 0) / length(r1_flags), 1)
    pct_rev  <- round(100 * sum(bitwAnd(r1_flags, 16L) >  0) / length(r1_flags), 1)
    message(sprintf("  Read1 forward : %.1f%%  |  Read1 reverse : %.1f%%",
                    pct_fwd, pct_rev))
    result <- if (pct_fwd >= 70) "forward stranded" else
      if (pct_rev >= 70) "reverse stranded" else "unstranded"
  }

  code <- switch(result,
                 "forward stranded" = 1, "reverse stranded" = 2, "unstranded" = 0)
  message(sprintf("\n  >> %s  (strandSpecific = %d in featureCounts/DESeq2)",
                  toupper(result), code))
  message("\nStep 3 complete.")
  return(invisible(result))
}

run_step4_split_strands <- function(bam_paths, paired_end = FALSE) {
  message("\n", strrep("=", 60))
  message("STEP 4: Splitting BAMs by strand")
  message(strrep("=", 60))

  result_list <- vector("list", length(bam_paths))

  for (i in seq_along(bam_paths)) {
    bam     <- bam_paths[i]
    prefix  <- sub("\\.bam$", "", bam)
    fwd_bam <- paste0(prefix, "_fwd.bam")
    rev_bam <- paste0(prefix, "_rev.bam")

    message(sprintf("\n  [%d/%d] %s", i, length(bam_paths), basename(bam)))

    if (!file.exists(paste0(bam, ".bai"))) {
      message("  Sorting & indexing: ", basename(bam))
      sorted_prefix <- sub("\\.bam$", "_sorted", bam)
      Rsamtools::sortBam(bam, sorted_prefix)
      file.rename(paste0(sorted_prefix, ".bam"), bam)
      Rsamtools::indexBam(bam)
    }

    if (!paired_end) {
      Rsamtools::filterBam(bam, fwd_bam,
                           param = Rsamtools::ScanBamParam(
                             flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                                                           isMinusStrand   = FALSE)))
      Rsamtools::filterBam(bam, rev_bam,
                           param = Rsamtools::ScanBamParam(
                             flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                                                           isMinusStrand   = TRUE)))
    } else {
      raw      <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(
        flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE),
        what = c("qname", "flag")))[[1]]
      is_read1 <- bitwAnd(raw$flag, 64L) > 0
      is_minus <- bitwAnd(raw$flag, 16L) > 0

      fwd_qnames <- unique(raw$qname[ is_read1 & !is_minus])
      rev_qnames <- unique(raw$qname[ is_read1 &  is_minus])

      Rsamtools::filterBam(bam, fwd_bam,
                           param  = Rsamtools::ScanBamParam(
                             flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE),
                             what = "qname"),
                           filter = S4Vectors::FilterRules(list(
                             fwd = function(x) x$qname %in% fwd_qnames)))

      Rsamtools::filterBam(bam, rev_bam,
                           param  = Rsamtools::ScanBamParam(
                             flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE),
                             what = "qname"),
                           filter = S4Vectors::FilterRules(list(
                             rev = function(x) x$qname %in% rev_qnames)))
    }

    Rsamtools::indexBam(fwd_bam)
    Rsamtools::indexBam(rev_bam)

    n_fwd <- length(Rsamtools::scanBam(fwd_bam,
                                       param = Rsamtools::ScanBamParam(what = "flag"))[[1]]$flag)
    n_rev <- length(Rsamtools::scanBam(rev_bam,
                                       param = Rsamtools::ScanBamParam(what = "flag"))[[1]]$flag)
    message(sprintf("  + strand : %d reads -> %s", n_fwd, basename(fwd_bam)))
    message(sprintf("  - strand : %d reads -> %s", n_rev, basename(rev_bam)))

    if (n_fwd == 0) warning("_fwd.bam is empty - _fwd.bedgraph will be blank. ",
                            "Run diagnose_bedgraph('", bam, "')")
    if (n_rev == 0) warning("_rev.bam is empty - _rev.bedgraph will be blank. ",
                            "Run diagnose_bedgraph('", bam, "')")

    result_list[[i]] <- list(fwd = fwd_bam, rev = rev_bam)
  }

  message("\nStep 4 complete.")
  return(invisible(result_list))
}

run_step5_bedgraph <- function(bam_paths, output_dir, paired_end = FALSE) {
  message("\n", strrep("=", 60))
  message("STEP 5: Converting BAM -> bedGraph")
  message(strrep("=", 60))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (i in seq_along(bam_paths)) {
    bam     <- bam_paths[i]
    bg_file <- file.path(output_dir, sub("\\.bam$", ".bedgraph", basename(bam)))
    message(sprintf("\n  [%d/%d] %s", i, length(bam_paths), basename(bam)))

    gal <- if (paired_end)
      GenomicAlignments::readGAlignmentPairs(bam, strandMode = 1)
    else
      GenomicAlignments::readGAlignments(bam)

    if (length(gal) == 0) {
      warning("No reads loaded from: ", basename(bam),
              " - bedGraph will be empty. Run diagnose_bedgraph('", bam, "')")
      file.create(bg_file)
      next
    }

    gr                  <- as(GenomicRanges::coverage(gal), "GRanges")
    colnames(mcols(gr)) <- "score"
    gr                  <- gr[gr$score > 0]

    if (length(gr) == 0) {
      warning("All coverage = 0 for: ", basename(bam),
              " - Run diagnose_bedgraph('", bam, "')")
      file.create(bg_file)
      next
    }

    message(sprintf("  Non-zero positions : %d", length(gr)))
    message(sprintf("  Coverage range     : %d - %d  (mean: %.1f)",
                    min(gr$score), max(gr$score), mean(gr$score)))
    rtracklayer::export(gr, bg_file, format = "bedGraph")
    message("  Saved: ", basename(bg_file))
  }
  message("\nStep 5 complete.")
}

# ---- Pipeline execution -------------------------------------------------

bam_paths <- run_step1_align(FASTQ_DIR, OUTPUT_DIR, ALIGN_INDEX, BAM_SUFFIX,
                             SAMPLE_NAMES, PAIRED_END, SEQ_TYPE)

run_step2_sort_index(bam_paths)

STRAND_RESULT <- run_step3_strandedness(bam_paths[1], paired_end = PAIRED_END)

if (SPLIT_STRAND) {
  strand_bam_list <- run_step4_split_strands(bam_paths, paired_end = PAIRED_END)
  fwd_bams <- sapply(strand_bam_list, `[[`, "fwd")
  rev_bams <- sapply(strand_bam_list, `[[`, "rev")
  run_step5_bedgraph(fwd_bams, OUTPUT_DIR, paired_end = PAIRED_END)
  run_step5_bedgraph(rev_bams, OUTPUT_DIR, paired_end = PAIRED_END)
} else {
  run_step5_bedgraph(bam_paths, OUTPUT_DIR, paired_end = PAIRED_END)
}

message("\n", strrep("=", 60))
message("PIPELINE COMPLETE")
message(strrep("=", 60))
message("  Index        : ", ALIGN_INDEX)
message("  Samples      : ", length(bam_paths))
message("  Strandedness : ", STRAND_RESULT)
message("  Strand split : ", SPLIT_STRAND)
message("  Output       : ", OUTPUT_DIR)
message("\n  Output files per sample:")
if (SPLIT_STRAND) {
  message("    <sample>_aligned.bam")
  message("    <sample>_aligned_fwd.bam  +  _fwd.bam.bai  (+ strand reads)")
  message("    <sample>_aligned_rev.bam  +  _rev.bam.bai  (- strand reads)")
  message("    <sample>_aligned_fwd.bedgraph              (+ strand coverage)")
  message("    <sample>_aligned_rev.bedgraph              (- strand coverage)")
} else {
  message("    <sample>_aligned.bam")
  message("    <sample>_aligned.bedgraph  (combined coverage)")
}
message(strrep("=", 60))
message("\nIf IGV shows nothing, run:")
message("  diagnose_bedgraph('path/to/sample.bam')")
message("  diagnose_bedgraph('path/to/sample.bedgraph')")
