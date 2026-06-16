#!/usr/bin/env Rscript
library(seqinr)
library(parallel)

#base_dir <- "/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/TES_analyses_Fig5/RNAfold/output/"
rnafold_path <- "/home/yip/viennaRNA/usr/bin/RNAfold"
num_cores <- 4 

#setwd(base_dir)

# Get the list of 12 gzipped fasta files
fasta_paths <- list.files(pattern = "fasta.gz")
file_names <- gsub("_up200_down200.fasta.gz", "", fasta_paths)

# --- HELPER FUNCTION ---
get_unpaired_probs <- function(ps_file, seq_len) {
  if (!file.exists(ps_file)) return(rep(NA, seq_len))
  
  ps_data <- readLines(ps_file)
  prob_lines <- ps_data[grep("^[0-9].*ubox$", ps_data)]
  
  p_paired <- rep(0, seq_len)
  
  if(length(prob_lines) > 0) {
    list_data <- strsplit(trimws(prob_lines), "\\s+")
    df <- as.data.frame(do.call(rbind, list_data), stringsAsFactors = FALSE)
    
    i_vals <- as.numeric(df[,1])
    j_vals <- as.numeric(df[,2])
    probs  <- as.numeric(df[,3])^2 
    
    for(k in 1:length(i_vals)) {
      if(i_vals[k] <= seq_len) p_paired[i_vals[k]] <- p_paired[i_vals[k]] + probs[k]
      if(j_vals[k] <= seq_len) p_paired[j_vals[k]] <- p_paired[j_vals[k]] + probs[k]
    }
  }
  return(p_paired) 
}

# --- WRAPPER FUNCTION FOR PARALLELIZATION ---
process_fasta_file <- function(idx) {
  current_fasta <- fasta_paths[idx]
  current_name  <- file_names[idx]
  
  tmp_dir <- paste0("tmp_", current_name)
  if(!dir.exists(tmp_dir)) dir.create(tmp_dir)
  
  fasta_data <- read.fasta(current_fasta, as.string = TRUE)
  
  all_probs <- list()      
  all_structures <- list()   
  
  # Move into temp directory to run RNAfold
  setwd(tmp_dir)
  
  for (i in 1:length(fasta_data)) {
    seq_id <- names(fasta_data)[i]
    sequence <- getSequence(fasta_data[[i]], as.string = TRUE)
    seq_len <- nchar(sequence)
    
    cmd <- sprintf("echo '>%s\n%s' | %s -p --noPS", seq_id, sequence, rnafold_path)
    rna_out <- system(cmd, intern = TRUE)
    
    mfe_struct <- if(length(rna_out) >= 3) rna_out[3] else NA
    centroid_struct <- if(length(rna_out) >= 6) rna_out[length(rna_out)] else NA
    
    all_structures[[seq_id]] <- data.frame(
      id = seq_id,
      mfe_structure = mfe_struct,
      centroid_structure = centroid_struct,
      stringsAsFactors = FALSE
    )
    
    clean_id <- strsplit(seq_id, " ")[[1]][1]
    ps_file <- paste0(clean_id, "_dp.ps")
    
    all_probs[[seq_id]] <- get_unpaired_probs(ps_file, seq_len)
    
    if(file.exists(ps_file)) file.remove(ps_file)
  }
  
  setwd("..") # Move back to output directory
  
  structure_df <- do.call(rbind, all_structures)
  write.table(structure_df, gzfile(paste0(current_name,"_structure_strings.tsv.gz")), 
              sep="\t", quote=FALSE, row.names=FALSE)
  
  center <- 201
  valley_idx <- (center - 5):(center + 14)
  body_idx   <- (center - 100):(center - 6)
  
  summary_scores <- data.frame(
    id = names(all_probs),
    mean_body_Pun = sapply(all_probs, function(x) mean(x[body_idx], na.rm=TRUE)),
    mean_valley_Pun = sapply(all_probs, function(x) mean(x[valley_idx], na.rm=TRUE))
  )
  
  summary_scores$depletion_ratio <- summary_scores$mean_valley_Pun / summary_scores$mean_body_Pun
  write.table(summary_scores, gzfile(paste0(current_name,"_structural_scoring.tsv.gz")), 
              sep="\t", quote=FALSE, row.names=FALSE)
  
  final_matrix <- do.call(rbind, all_probs)
  saveRDS(final_matrix, paste0(current_name,"_full_prob_matrix.rds"))
  
  # Cleanup temp directory
  unlink(tmp_dir, recursive = TRUE)
  
  return(paste0("Success: ", current_name))
}

# ---  EXECUTION ---
cat(paste0("Starting Parallel Processing on ", num_cores, " cores...\n"))

mclapply(1:length(fasta_paths), process_fasta_file, mc.cores = num_cores)

cat("\nAll 12 files processed successfully.\n")

