reorderIdents <- function(object, ident_type, desired_order, built_in = FALSE){

  Seurat::Idents(object) <- object[[ident_type]]
  z <- object@active.ident; levels(z)

  if (built_in){
    if (ident_type=="stim"){
      z_desired_order <- c("Ctrl","LPS 500ug 3h","LPS 500ug 6h","LPS 500ug 1d","LPS 500ug 2d","LPS >3w","LPS 10mg 3h","LPS 10mg 6h",
                           "LPS 50mg 3h","LPS 50mg 6h","Poly(i:c) 10mg 3h","Poly(i:c) 10mg 6h","Poly(i:c) 20mg 3h","Poly(i:c) 20mg 6h",
                           "Poly(i:c) >3w","TNF-α 500ug 6h")
    } else if (ident_type=="cell_type"){
      z_desired_order <- c("Somatotropes", "Corticotropes", "Lactotropes", "Gonadotropes", "Melanotropes", "Thyrotropes", "Pituicytes", "Pericytes",
                           "Endothelial cells", "Stem cells", "Pou1f1 progenitors", "Folliculostellate cells", "White blood cells", "Red blood cells")
    } else if (ident_type=="treat"){
      z_desired_order <- c("Saline","LPS","Poly(i:c)","TNF-α")
    }
  } else {
    z_desired_order <- desired_order
  }

  z <- factor(z,levels(z)[ match(z_desired_order,levels(z)) ])  # Change factor level orders
  object[[ident_type]] <- z
  return(object)

}
