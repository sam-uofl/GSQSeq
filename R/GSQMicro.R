#' This function perofrms gene set analysis with Quantitative Trait Loci for Gene expreesion data derrived from Microarray studies through integrating differential expression (DE) scores of genes.
#'
#' @param ExpData N by M matrix of expression values, rows as genes (N: total number of genes) and columns as samples (M: Total samples).
#' @param class Vector of 1 and 0 of length M shows the class labels of samples (1: case and 0: control).
#' @param geneLoc N by 3 dataframe/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene set: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs in their respective chromosomes.
#' @param size size of the selected gene set to be analyzed with the QTL, e.g. 100, 500, ...
#' @param QTLData Q by 3 dataframe/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs in their respective chromosomes.
#' @param method A character representing the Differential Expression (DE) analysis method to be used for DE analysis of expression data. It must be either T-test or F-scores or MRMR or SVM.
#'
#' @return This functions rturns the statistical QTL enrichment significance value of the gene set of size 'size' along with test statistic.
#'
#' @author Samarendra Das
#'
#'
#' @importFrom stats pnorm
#' @importFrom stats t.test
#' @importFrom stats var
#' @importFrom BootMRMR Weights.mrmr
#' @importFrom e1071 svm
#' @importFrom GSAQ totqtlhit
#'
#' @export

GSQMicro <- function (ExpData, class, geneLoc, size, QTLData, method){

  if(!is.matrix(ExpData) & !is.data.frame(ExpData) & class(ExpData)[1] != "dgCMatrix")
    stop("Wrong input data type of 'Exp Data'")
  if(sum(is.na(ExpData)) > 0)
    stop("NAs detected in input 'Expression Data'");gc();
  if(sum(ExpData < 0) > 0)
    stop("Negative values detected in input 'Exp Data'");gc();
  if(all(ExpData == 0))
    stop("All elements of input 'Exp Data' are zeros");gc();
  if(any(colSums(ExpData) == 0))
    warning("Library size of zeros detected in 'Exp Data'");gc();
  #if(!is.factor(group) )
  #stop("Data type of 'group' is not factor")
  if(max(class) != 1)
    stop("Levels number of 'class' is not two")
  if(table(class)[1] < 2 | table(class)[2] < 2)
    stop("Too few samples (< 2) in a group")
  if(ncol(ExpData) != length(class))
    stop("Length of 'class' must equal to column number of 'Exp Data'")

  if(!is.character(method))
    stop("Data type of 'method' is not character")
  if(!is.element(method, c("Ttest", "F", "MRMR", "SVM")))
    stop("method mlust be either from  'Ttest', 'F-score', 'MRMR', 'SVM'")

  if(is.null(method)) method <- "Ttest"

  #ExpData = Expdata
  #QTLData = QTLdata
  N <- nrow(ExpData)               ###number of genes####
  M <- ncol(ExpData)          ###number of samples###
  class <- as.vector(class)
  if (method == "Ttest")
  {
    # ranking <- vector(length=n)
    idx <- which(class == 1)  # indexing of positive samples
    idy <- which(class == 0) # indexing of negative samples
    #print(idx)
    B=vector(mode="numeric", N)
    for(i in 1:N){
      t.mes <-t.test(ExpData[i, idx], ExpData[i, idy])$statistic  #####F-Score
      B[i] <- t.mes
      names(B) <- rownames(ExpData)
    }

    B = sort(abs(B), decreasing = TRUE, index.return = FALSE)
  }

  #########method
  if(method == "F"){
    B = vector(mode = "numeric", N)
    idx <- which(class == 1)  # indexing of positive samples
    idy <- which(class == 0) # indexing of negative samples
    g <- as.matrix(ExpData)
    for (i in 1:N) {
      f.mes <- (((mean(g[i, idx]) - mean(g[i, ]))^2) + ((mean(g[i,
                                                                idy]) - mean(g[i, ]))^2))/(var(g[i, idx]) + var(g[i,
                                                                                                                  idy]))
      B[i] <- f.mes
    }
    names(B) <- row.names(ExpData)
    B = sort(abs(B), decreasing = TRUE, index.return = FALSE)

  }

  if(method == "MRMR"){
    class <- ifelse(class == 0, -1, 1)
    B <- Weights.mrmr(x = ExpData, y = class)
    B <- as.vector(B)
    names(B) <- row.names(ExpData)
    B <- sort(abs(B), decreasing = TRUE, index.return = FALSE)
  }

  if(method == "SVM"){
    if (!requireNamespace("e1071", quietly = TRUE)) {
      stop("Package \"e1071\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    y <- as.factor(class)
    genes <- rownames(ExpData)
    x1 <- as.matrix(ExpData)
    #usethis::use_package("e1071")
    SvmModel <- svm(t(x1), y, cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" )
    a <- SvmModel$coefs
    b <- SvmModel$SV
    B <- abs((t(a)%*%b)) # absolute value of SVM weights
    names(B) <- genes
    B = sort(abs(B), decreasing = TRUE, index.return = FALSE)
  }

  #size <- 500
  GeneSel <- B[1:size]
  #length( GeneSel)
  GenSelpr <- B[-(1:size)]
  #length(GenSelpr)
  #GenSelpr[1:10]
  gene.sel <- names(GeneSel)
  geneqtl <- GeneQTL (geneset = gene.sel, genelist = geneLoc, qtl = QTLData)
  geneqtl <- geneqtl$Mat1[, c(1:5)]
  #dim(geneqtl)
  idd <- GeneSel[match(as.vector(geneqtl[,1]), names(GeneSel))]
  #sum(idd)
  gene.sel.pr <- names(GenSelpr)
  geneprqtl <- GeneQTL (geneset = gene.sel.pr, genelist = geneLoc, qtl = QTLData)
  #dim(geneprqtl)
  geneprqtl <- geneprqtl$Mat1[, c(1:5)]
  iddpr <- GenSelpr[match(as.vector(geneprqtl[,1]), names(GenSelpr))]
  #length(iddpr)
  #sum(iddpr)
  Diff <- sum(idd) - sum(iddpr)

  ###expected value of diff scores of genes in gene set
  exp.genset <- mean(as.numeric(GeneSel))
  var.genset <- var(as.numeric(GeneSel))
  exp.genset2 <- var.genset + exp.genset^2
  ########expected value and variance of qtlhit genes in geneset
  V <- totqtlhit(genelist = geneLoc, qtl = QTLData)  ######total hits in gene space
  V <- as.numeric(V[1])
  v <- nrow(geneqtl)
  N <- nrow(ExpData)
  n <- length(GeneSel)
  exp.ngq <- as.numeric(n * V / N)
  var.ngq <- n * V *(N - V) * (N - n) / ((N - 1) * N^2)

  ###########expected value and variance of diff
  exp.diff <- 2 * exp.genset *  exp.ngq - v * exp.genset
  var.diff <- 4 * var.genset/(v-1) * exp.ngq * (v - exp.ngq) - var.ngq + exp.genset2 * var.ngq
  var.diff <- abs(var.diff)
  ######distribution of disff stat
  Z <- (Diff - exp.diff)/ (var.diff)^0.5
  pval <- pnorm(abs(Z), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  out <- list(QTLHitGeneSet = idd, statistic = Diff, p.value = pval)
  class(out) <- c("QTLhit genes in geneset", "Test  Statistic", "p-value")
  return(out)
}
