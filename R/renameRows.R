#' @title Rename Row Names (Ensembl ID to Symbol)
#'
#' @description For generic analysis of data not generated commercial pipelines (i.e. SMART-SEQ2), we are usually given a raw count matrix,
#' where rows are Ensembl IDs. To rename the rows from Ensembl IDs to gene Symbols, a few things need to be considered:
#'
#' First, row names of Seurat object can NOT be modified, so the renaming operation MUST be done on original expression matrix, usually loaded
#' as a dataframe.
#' Second, different Ensembl IDs could map to the same gene Symbol, causing duplicated row names, which is forbidden in dataframe. To
#' circumvent this, we attach a digit to the Symbols.
#'
#' @param data A dataframe, rows are features named as Ensembl IDs, columns are cells.
#' @param anno Annotation table. Used to search mapping relationships between Ensembl IDs and Symbols.
#' @return A dataframe, rows are features renamed as gene Symbols.
#' @export
#'
renameRows <- function(df, anno){
  df <- dplyr::slice(df, -which(!rownames(df) %in% anno$ensembl_id))  # df may contain rows not in the annotation table. Remove these rows.
  rn <- rownames(df)
  rn <- anno[match(rn,anno$ensembl_id),"symbol"]
  dup <- rn[duplicated(rn)]
  print(paste0("The following genes are duplicated: ", paste0(dup, collapse = ", "), "."))
  for (d in dup){
    count = 0
    for (i in 1:length(rn)){
      if (rn[i] == d) {
        count = count+1
        if (count>1){
          rn[i] = paste0(rn[i],"-",count)
        }
      }
    }
  }
  rownames(df) <- rn
  return(df)
}
