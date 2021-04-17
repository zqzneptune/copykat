#' annotate genes with reference to hg38.
#'
#' @param mat data matrix; genes in rows; cell names in columns.
#' @param ID.type gene id type: Symbol or Ensemble.
#' @param full.anno.hg20 annotation file for all known genes, automatically loaded in copycat.
#'
#' @return annotations of each genes in rows with chrom and positions.
#'
#' @examples
#' test.anno.mat <- annotateGenes.hg20(mat=matx, ID.type="ENSEMBLE_id", full.anno.hg20 = full.anno.hg20)
#' @export
annotateGenes.hg20 <- function(mat, ID.type="S"){
  print("start annotation ...")

  if(substring(ID.type,1,1) %in% c("E", "e")){
    shar <- intersect(rownames(mat), full.anno.hg20$ensembl_gene_id)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno.hg20[which(as.vector(full.anno.hg20$ensembl_gene_id) %in% shar),]
    anno <- anno[!duplicated(anno$gene_symbol),]
    anno <- anno[order(match(anno$ensembl_gene_id, rownames(mat))),]
    data <- cbind(anno, mat)

  }else if(substring(ID.type,1,1) %in% c("S", "s")) {

    shar <- intersect(rownames(mat), full.anno.hg20$gene_symbol)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno.hg20[which(as.vector(full.anno.hg20$gene_symbol) %in% shar),]
    anno <- anno[!duplicated(anno$gene_symbol),]
    anno <- anno[order(match(anno$gene_symbol, rownames(mat))),]
    data <- cbind(anno, mat)
  }
}

