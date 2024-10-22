% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iCNplot.R
\name{iCNplot}
\alias{iCNplot}
\title{Plot post-segmentation copy number profiles of integer values}
\usage{
iCNplot(iCNmat, ref, Gini, annotation = NULL,
        plot.dendrogram = TRUE, filename)
}
\arguments{
\item{iCNmat}{inferred integer copy-number matrix by SCOPE,
with each column being a cell and each row being a genomic bin}

\item{ref}{GRanges object after quality control procedure}

\item{Gini}{vector of Gini coefficients for each cell,
with the same order as that of cells in columns of \code{iCNmat}}

\item{annotation}{vector of annotation for each cell,
with the same order as that of cells in columns of \code{iCNmat}.
Default is \code{NULL}.}

\item{plot.dendrogram}{logical, whether to plot the dendrogram.
Default is \code{TRUE}.}

\item{filename}{name of the output png file}
}
\value{
png file with integer copy-number profiles across single cells
with specified annotations
}
\description{
Show heatmap of inferred integer copy-number profiles
by SCOPE with cells clustered by hierarchical clustering
}
\examples{
Gini = getGini(Y_sim)
iCNplot(iCNmat = iCN_sim,
        ref = ref_sim,
        Gini = Gini,
        filename = 'iCNplotdemo')

}
\author{
Rujin Wang \email{rujin@email.unc.edu}
}
