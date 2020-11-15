#'
#' This function enables to find the list of the selected genes overlapped with QTL regions.
#'
#' @param geneset vector of characters representing the names of genes/ gene ids selected from the whole gene list/space by using a gene selection method.
#' @param genelist N by 3 dataframe/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene set: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs in their respective chromosomes.
#' @param qtl Q by 3 dataframe/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs in their respective chromosomes.
#'
#' @description The function provides the list of the selected genes overlapped Quantitative Trait Loci (QTL) regions along with their genomic positions.
#' @retun This functio returns the list of the QTL-hit genes in the selected gene set.
#'
#' @author Samarendra Das
#'
#' @export


GeneQTL <- function (geneset, genelist, qtl)
{
  this.call = match.call()
  if ((!class(geneset) == "character")) {
    warning("gene set must be a vector of gene names")
  }
  if (!class(genelist) == "data.frame") {
    warning("genelist must be a dataframe with row names as gene ids and three columns as chromosome number, start and stop positions")
  }
  if (length(geneset) > nrow(genelist)) {
    warning("geneset must be a sub set of genelist")
  }
  genenames <- rownames(genelist)
  temp <- match(unlist(geneset), genenames)
  if (sum(is.na(temp)) == length(temp)) {
    stop("IDs in the gene set and genelist are not matching. Ensure same type of IDs in both the sets")
  }
  if (sum(!is.na(temp))/length(temp) < 0.05) {
    stop("Fewer than 5% of genes in the genesets appear in the dataset. Make sure\that gene identifiers in dataset are Gene symbols")
  }
  qtl.start <- as.vector(qtl[, 2])
  qtl.stop <- as.vector(qtl[, 3])
  m <- nrow(qtl)
  n <- nrow(geneset)
  geneset <- as.vector(geneset)
  genenames <- rownames(genelist)
  temp1 <- match(unlist(geneset), genenames)
  tempp <- temp1[!is.na(temp1)]
  select.gen <- genelist[tempp, ]
  chr.gen <- select.gen[, 1]
  chr.qtl <- qtl[, 1]
  gene.start <- as.vector(select.gen[, 2])
  gene.stop <- as.vector(select.gen[, 3])
  chr.slt <- as.numeric(names(table(chr.qtl)))
  selgenenam <- rownames(select.gen)
  Mat <- NULL
  Mat1 <- NULL
  for (i in 1:length(selgenenam)) {
    id <- as.vector(which(chr.qtl == chr.gen[i]))
    if (length(id) > 0) {
      a <- matrix("NULL", 1, length(id))
      b <- matrix("NULL", 1, length(id))
      for (j in 1:length(id)) {
        if (gene.start[i] >= qtl.start[id[j]] & gene.stop[i] <=
            qtl.stop[id[j]]) {
          a <- cbind(rownames(select.gen[i, ]), rownames(qtl[id[j],
                                                             ]))
          b <- cbind(rownames(select.gen[i, ]), select.gen[i,
                                                           ], rownames(qtl[id[j], ]), qtl[id[j], ])
          Mat <- rbind(Mat, a)
          Mat1 <- rbind(Mat1, b)
        }
      }
    }
  }
  colnames(Mat) <- c("Gene ID", "QTL ID")
  Mat <- as.data.frame(Mat)
  rownames(Mat1) <- NULL
  colnames(Mat1) <- c("Gene ID", "Chr. No.", "Start (bp)",
                      "End (bp)", "QTL ID", "Chr. No", "Start (bp)",
                      "End (bp)")
  res <- list(Mat = Mat, Mat1 = Mat1)
  class(res) <- c("List of qtlhit genes", "Genomic positions of qtlhit genes")
  #res <- Mat1[, 1:5]
  return(res)
}

###########example
#GeneQTL(geneset=geneset, genelist=geneLoc, qtl=QTLData)
