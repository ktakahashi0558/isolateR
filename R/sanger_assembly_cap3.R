# -------------------------------------------------------------------------
# 内部ヘルパー関数 (Internal Helper Function for CAP3)
# -------------------------------------------------------------------------
assemble_with_cap3 <- function(dna_string_set,
                               cap3_path = "cap3",
                               percent_identity = 90,
                               overlap_length = 20,
                               similarity_cutoff = 300,
                               max_qscore_sum = 5000) {
  
  # 一時ファイルのベース名を生成
  temp_base <- tempfile()
  input_fasta <- paste0(temp_base, ".fasta")
  
  # リードをFASTAファイルに書き出す
  Biostrings::writeXStringSet(dna_string_set, filepath = input_fasta)
  
  # CAP3コマンドを構築
  cap3_command <- paste(
    cap3_path,
    shQuote(input_fasta),
    "-p", percent_identity,
    "-o", overlap_length,
    "-s", similarity_cutoff,
    "-d", max_qscore_sum
  )
  
  cat("Running CAP3...\n")
  cat("Command:", cap3_command, "\n")
  
  # CAP3は大量の出力ファイルを生成するため、作業ディレクトリを一時的に変更
  original_wd <- getwd()
  temp_dir <- dirname(temp_base)
  setwd(temp_dir)
  
  # CAP3を実行
  system(cap3_command, intern = TRUE, ignore.stderr = TRUE)
  
  # 作業ディレクトリを元に戻す
  setwd(original_wd)
  
  # 結果ファイル（コンティグ）のパス
  contigs_file <- paste0(input_fasta, ".cap.contigs")
  
  final_contig <- Biostrings::DNAString("")
  
  if (file.exists(contigs_file) && file.info(contigs_file)$size > 0) {
    contigs <- Biostrings::readDNAStringSet(contigs_file)
    if (length(contigs) > 0) {
      # 最長のコンティグを最終結果とする
      final_contig <- contigs[[which.max(width(contigs))]]
    }
  }
  
  # CAP3が生成した全一時ファイルを削除
  unlink(list.files(path = temp_dir, pattern = basename(temp_base), full.names = TRUE))
  
  return(final_contig)
}


# -------------------------------------------------------------------------
# メイン関数 (Exported Main Function)
# -------------------------------------------------------------------------
#' @title Assemble Sanger sequencing reads robustly using CAP3
#' @description This function takes a data frame of Sanger reads, groups them by a common 
#' filename prefix, and assembles each group into a consensus sequence using the 
#' external command-line tool CAP3.
#'
#' @param input A character string specifying the path to the input CSV file. The CSV must
#' contain a 'filename' column and a 'seqs_trim' column with the DNA sequences.
#' @param suffix A regular expression string used to remove suffixes from filenames to identify groups.
#' @param percent_identity A numeric value for the minimum percent identity for overlaps in CAP3 (default: 95).
#' @param overlap_length A numeric value for the minimum overlap length in CAP3 (default: 40).
#' @param cap3_path A character string specifying the path to the CAP3 executable.
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A data frame or an isoQC object containing the assembly results.
#'
#' @importFrom Biostrings DNAString DNAStringSet readDNAStringSet width writeXStringSet
#' @importFrom utils read.csv write.csv
#' @importFrom stringr str_count
#'
#' @export
sanger_assembly_cap3 <- function(input = NULL,
                                  suffix = "_F.ab1|_R.ab1",
                                  percent_identity = 95,
                                  overlap_length = 20,
                                  cap3_path = "cap3",
                                  verbose = TRUE) {
  
  if (!requireNamespace("Biostrings", quietly = TRUE)) stop("Biostrings required")
  if (is.null(input)) stop("input must be provided")
  
  input.df <- utils::read.csv(input, stringsAsFactors = FALSE)
  input.df$no_suffix <- sub(paste0("(", suffix, ")$"), "", input.df$filename)
  input.df$row_counter <- seq_len(nrow(input.df))
  
  counts <- table(input.df$no_suffix)
  paired_names <- names(counts[counts >= 2])
  unpaired_names <- names(counts[counts < 2])
  
  paired.seqs <- input.df[input.df$no_suffix %in% paired_names, , drop = FALSE]
  unpaired.seqs <- input.df[input.df$no_suffix %in% unpaired_names, , drop = FALSE]
  
  collector.list <- vector("list", length(paired_names))
  
  for (ii in seq_along(paired_names)) {
    grp <- paired_names[ii]
    paired.seqs.filt <- paired.seqs[paired.seqs$no_suffix == grp, , drop = FALSE]
    
    original_seqs_char <- as.character(paired.seqs.filt$seqs_trim)
    sanitized_seqs_char <- gsub("[^ACGTNMRWSYKVHDB-]", "N", toupper(original_seqs_char), perl = TRUE)
    all_reads_in_group <- Biostrings::DNAStringSet(sanitized_seqs_char)
    
    final_contig <- Biostrings::DNAString("")
    failed.aln <- TRUE
    
    # S4オブジェクトには width() を使う
    if (length(all_reads_in_group) < 2 || all(width(all_reads_in_group) == 0)) {
      message("Group ", grp, " has fewer than 2 reads or only empty sequences. Skipping assembly.")
      out_row <- paired.seqs.filt[1, , drop = FALSE]
      out_row$decision <- "Fail"
    } else {
      tryCatch({
        final_contig <- assemble_with_cap3(all_reads_in_group, cap3_path, percent_identity, overlap_length)
        failed.aln <- (length(final_contig) == 0)
        if (!failed.aln) {
          message("Group ", grp, " successfully assembled via CAP3. Contig length: ", length(final_contig))
        } else {
          message("Group ", grp, " assembly failed with CAP3. No contig was formed.")
        }
      }, error = function(e) {
        message("Group ", grp, " failed during CAP3 execution. Error: ", e$message)
        final_contig <<- Biostrings::DNAString("")
        failed.aln <<- TRUE
      })
    }
    
    out_row <- paired.seqs.filt[1, , drop = FALSE]
    out_row$filename <- paste0(grp, "_consensus")
    out_row$decision <- ifelse(failed.aln, "Fail", "Pass")
    out_row$seqs_trim <- as.character(final_contig)
    out_row$length_trim <- length(final_contig)
    out_row$Ns_trim <- stringr::str_count(out_row$seqs_trim, "N")
    out_row$aln <- ifelse(failed.aln, "fail", "pass")
    
    collector.list[[ii]] <- out_row
  }
  
  collector.df <- do.call(rbind, collector.list)
  
  if (nrow(unpaired.seqs) > 0) {
    missing_cols <- setdiff(names(collector.df), names(unpaired.seqs))
    for (col in missing_cols) {
      unpaired.seqs[[col]] <- NA
    }
    unpaired.seqs <- unpaired.seqs[, names(collector.df)]
  }
  
  final.df <- rbind(collector.df, unpaired.seqs)
  final.df <- final.df[order(final.df$row_counter), ]
  
  csv.out.path <- paste0(sub("\\.csv$", "", input), "_consensus.csv")
  utils::write.csv(final.df, csv.out.path, row.names = FALSE)
  
  return(final.df)
}
