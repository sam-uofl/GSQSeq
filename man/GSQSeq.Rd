% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSQSeq.R
\name{GSQSeq}
\alias{GSQSeq}
\title{This function performs analysis of selected genes from the RNA-seq data with QTL annotation data.}
\usage{
GSQSeq(CountData, class, geneLoc, size, QTLData, method)
}
\arguments{
\item{CountData}{N by M matrix of reads count data (N: Total number of genes; M: Total number of samples), rows represent the genes and columns represent samples/libraries.}

\item{class}{M by 1 vector of 0 and 1 representing the class labels of the samples, i.e. 1 for case and 0 for control.}

\item{geneLoc}{N by 3 data frame/ matrix (genes/gene ids as row names); where, N represents the number of genes in the whole gene set: first coloumn represnting the chromosomal location of genes: second coloumn representing the start position of genes in terms of basepairs: third coloumn representing the end position of genes in terms of basepairs in their respective chromosomes.}

\item{size}{size of the selected gene set to be analyzed with the QTL, e.g. 100, 500, ...}

\item{QTLData}{Q by 3 data frame/matrix (qtl names/qtl ids as row names);where, Q represents the number of qtls: first coloumn represnting the chromosomal location of qtls: second coloumn representing the start position of qtls in terms of basepairs: third coloumn representing the end position of qtls in terms of basepairs in their respective chromosomes.}

\item{method}{A character representing the Differential Expression (DE) analysis method to be used for DE analysis of RNA-seq data. It must be either DESeq or edgeR.}
}
\value{
This function returns the output from the gene set analysis of the selected genes with QTL data.
}
\description{
This function performs analysis of selected genes from the RNA-seq data with QTL annotation data.
}
\author{
Samarendra Das
}
