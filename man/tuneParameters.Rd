% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc.R
\name{tuneParameters}
\alias{tuneParameters}
\title{Optimize Resolution and Gene Number Parameters for Cell Type Annotation}
\usage{
tuneParameters(srt, resolution = numeric(), N = integer(), spc)
}
\arguments{
\item{srt}{Seurat object, the input data object to be analyzed.}

\item{resolution}{numeric vector, a vector of resolution values to be tested in \code{Seurat::FindClusters()}.}

\item{N}{integer vector, a vector of values indicating the number of top differential genes to be used for matching in \code{matchCellMarker2()}.}

\item{spc}{character, the species parameter for the \code{matchCellMarker2()} function, specifying the organism.}
}
\value{
A list of ggplot2 objects, each representing a UMAP plot generated with a different combination of resolution and N parameters.
}
\description{
This function tunes the \code{resolution} parameter in \code{Seurat::FindClusters()} and the number of top differential genes (\code{N}) to obtain different cell type annotation results. The function generates UMAP plots for each parameter combination, allowing for a comparison of how different settings affect the clustering and annotation.
}
