#'
#' This function performs analysis of selected genes from the RNA-seq data with QTL annotation data.
#' @param CountData N by M matrix of reads count data (N: Total number of genes; M: Total number of samples), rows represent the genes and columns represent samples/libraries.
#' @param class M by 1 vector of 0 and 1 representing the class labels of the samples, i.e. 1 for case and 0 for control.
#' @param geneLoc N by 3 data frame/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene set: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs in their respective chromosomes.
#' @param size size of the selected gene set to be analyzed with the QTL, e.g. 100, 500, ...
#' @param QTLData Q by 3 data frame/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs in their respective chromosomes.
#' @param method A character representing the Differential Expression (DE) analysis method to be used for DE analysis of RNA-seq data. It must be either DESeq or edgeR.
#'
#' @return This function returns the output from the gene set analysis of the selected genes with QTL data.
#'
#' @author Samarendra Das
#'
#' @importFrom stats pnorm
#' @importFrom stats var
#' @importFrom stats model.matrix
#' @import DESeq2
#' @import edgeR
#' @importFrom S4Vectors DataFrame
#'
#' @export
#'
GSQSeq <- function (CountData, class, geneLoc, size, QTLData, method){
  if(!is.matrix(CountData) & !is.data.frame(CountData) & class(CountData)[1] != "dgCMatrix")
    stop("Wrong input data type of 'Count Data'")
  if(sum(is.na(CountData)) > 0)
    stop("NAs detected in input 'Count Data'");gc();
  if(sum(CountData < 0) > 0)
    stop("Negative values detected in input 'Count Data'");gc();
  if(all(CountData == 0))
    stop("All elements of input 'Count Data' are zeros");gc();
  if(any(colSums(CountData) == 0))
    warning("Library size of zeros detected in 'Counts Data'");gc();
  #if(!is.factor(group) )
  #stop("Data type of 'group' is not factor")
  if(max(class) != 1)
    stop("Levels number of 'class' is not two")
  if(table(class)[1] < 2 | table(class)[2] < 2)
    stop("Too few samples (< 2) in a group")
  if(ncol(CountData) != length(class))
    stop("Length of 'class' must equal to column number of 'Count Data'")

  if(!is.character(method))
    stop("Data type of 'method' is not character")
  if(!is.element(method, c("DEseq", "edgeR")))
    stop("method mlust be either from  'DEseq', 'edgeR'")

  if(is.null(method)) method <- "edgeR"

  if(method == "DESeq"){
    group <- as.factor(class)
    CountData <- as.matrix(CountData)
    #require(DESeq2)
    dds <- DESeqDataSetFromMatrix(CountData, DataFrame(group), ~ group)
    ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
    resLRT <- results(ddsLRT)
    B <- resLRT$stat
    names(B) <- rownames(resLRT)
    remove(resLRT)
  }
  #countData <- countData[1:1000, ]
  if (method == "edgeR"){
    group <- as.factor(class)
    y <- DGEList(counts = CountData, group = group)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, model.matrix(~group))
    fit <- glmFit(y,design = model.matrix(~group))
    lrt <- glmLRT(fit, coef=2:length(levels(group)))
    lrt <- lrt$table
    B <- lrt[, 3]
    names(B) <- rownames(lrt)
    remove(lrt)
  }

B = sort(abs(B), decreasing = TRUE, index.return = FALSE)
#size <- 500
GeneSel <- B[1:size]
#GeneSel[1:10]
#length( GeneSel)
GenSelpr <- B[-(1:size)]
#length(GenSelpr)
#GenSelpr[1:10]
gene.sel <- names(GeneSel)
geneqtl <- GeneQTL (geneset = gene.sel, genelist = geneLoc, qtl = QTLData)
#geneqtl <- geneqtl[, 1:5]
geneqtl <- geneqtl$Mat1[, c(1:5)]
dim(geneqtl)
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
N <- nrow(CountData)
M <- ncol(CountData)
n <- length(GeneSel)
exp.ngq <- as.numeric(n * V / N)
var.ngq <- n * V *(N - V) * (N - n) / ((N - 1) * N^2)

###########expected value and variance of diff
exp.diff <- 2 * exp.genset *  exp.ngq - v * exp.genset
var.diff <- 4 * var.genset/(v-1) * exp.ngq * (v - exp.ngq) - var.ngq + exp.genset2 * var.ngq
var.diff <- abs(var.diff)
######distribution of disff stat
Z <- (Diff - exp.diff)/ (var.diff)^0.5
pval <- pnorm(abs(Z), mean = 0, 1, lower.tail = FALSE, log.p = FALSE)
out <- list(QTLHitGeneSet = idd, statistic = Diff, p.value = pval)
class(out) <- c("QTLhit genes in geneset", "Test  Statistic", "p-value")
return(out)
}
