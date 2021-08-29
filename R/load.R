#' @title Load gencode annotation from file
#'
#' @description Load gene name conversion table, extracted from gencode vM21 annotation gtf.
#'
#' @export
#'
loadGencodeAnno <- function(strip_id = TRUE){
  table <- read.table("~/Documents/Work/scRNAseq/common_resources/gencode.vM21.annotation.tab",
                      sep = "\t", col.names = c('ensembl_id','symbol'))
  if(strip_id){
    table$stable_id <- sub(pattern = "\\.\\d+", "",x = table$stable_id)
  }
  return(table)
}

loadIUPHARDatabase <- function(){
  table <- read.csv("~/Documents/Work/scRNAseq/common_resources/targets_and_families.csv")
  return(table)
}

loadEnsemblAnno <- function(){
  table <- read.table("~/Documents/Work/scRNAseq/common_resources/MM.GRCm38.102.annotation.tab",
                      sep = "\t", col.names = c('ensembl_id','symbol'))
  return(table)
}
