% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scHolography.R
\name{dataAlign}
\alias{dataAlign}
\title{Data Integration}
\usage{
dataAlign(
  low.res.sp,
  high.res.sp,
  stProcessed = F,
  scProcessed = F,
  nFeatureToUse = 3000,
  whichReference = 2,
  nPCtoUse = 32,
  future.size = 4000
)
}
\arguments{
\item{low.res.sp}{Seurat Object of ST (low resolution) data}

\item{high.res.sp}{Seurat Object of SC (high resolution) data}

\item{stProcessed}{Is low.res.sp processed? Default is False}

\item{scProcessed}{Is low.res.sp processed? Default is False}

\item{nFeatureToUse}{Number of features to use for downstream analysis (e.g., PCA). Default is 3000}

\item{whichReference}{Use which dataset as the reference for Seurat integration. Default using SC data (whichReference=2). Set whichReference=1 to use ST data as the integration reference}

\item{nPCtoUse}{Number of PCs used for analysis. Default is 32}
}
\description{
Data Integration
}
