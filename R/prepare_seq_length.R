#' Calculate lengths of Transcripts' peptide sequences if exists
#'
#' @param fastaFileDir fasta file of ENST sequences
#' @importFrom seqinr read.fasta
#' @return dataframe
#' @export
#'
prepare_seq_length <- function(fastaFileDir) {

  biomart_protein_seq <- seqinr::read.fasta(file = fastaFileDir,
                                    seqtype = c("AA"), as.string = TRUE,
                                    bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong,
                                    endian = .Platform$endian, apply.mask = TRUE)

  df_biomart_seq <-  as.data.frame(do.call(rbind, biomart_protein_seq))
  head(df_biomart_seq, 2)

  genes_transcripts <- strsplit(rownames(df_biomart_seq), "\\|")
  genes <- sapply(genes_transcripts, `[`, 1)
  transcripts <- sapply(genes_transcripts, `[`, 2)

  # Add these as new columns
  df_biomart_seq$ENSG <- genes
  df_biomart_seq$ENST <- transcripts

  # Calculate the number of characters in the V1 column
  df_biomart_seq$Length <- nchar(df_biomart_seq$V1)
  df_biomart_seq_selection <- df_biomart_seq[, c('ENSG', 'ENST', 'Length')]

  return(df_biomart_seq_selection)
}
