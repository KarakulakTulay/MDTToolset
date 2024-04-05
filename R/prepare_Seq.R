#' Title
#'
#' @param fastaFileDir fasta file of ENST sequences
#' @param ENST_geneNamesDir ENST Gene Name/ENSG matching
#' @importFrom Biostrings readDNAStringSet
#' @return dataframe
#' @export
#'
prepare_Seq <- function(fastaFile, ENSG_ENST_ENSP) {

  ENST_sequences <- fastaFile

  ENST_sequences_df <- function(x) data.frame(width=Biostrings::width(x), seq=as.character(x), names=names(x))
  ENST <- gsub('\\..*', '', names(ENST_sequences))

  ENST_df <- ENST_sequences_df(ENST_sequences)
  ENSTs_sequences <- cbind(as.data.frame(ENST), ENST_df[,2])

  # read ENST HGNC
  ENSG_ENST_ENSP <- ENSG_ENST_ENSP
  ENSG_ENST <- ENSG_ENST_ENSP[,1:2]
  colnames(ENSG_ENST) <- c('ENSG', 'ENST')

  ## ensg_enst_seq
  ensg_enst_seq <- merge(ENSG_ENST, ENSTs_sequences, by = 'ENST')
  colnames(ensg_enst_seq) <-  c('ENST', 'ENSG', 'ENST_Seq')
  return(ensg_enst_seq)
}
