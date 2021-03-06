% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/familias2ped.R
\name{Familias2ped}
\alias{Familias2ped}
\alias{readFamiliasLoci}
\title{Convert \code{Familias} objects to \code{ped} format}
\usage{
Familias2ped(familiasped, datamatrix, loci, matchLoci = FALSE)

readFamiliasLoci(loci)
}
\arguments{
\item{familiasped}{A \code{FamiliasPedigree} object or a list of such.}

\item{datamatrix}{A data frame with two columns per marker (one for each
allele) and one row per individual.}

\item{loci}{A \code{FamiliasLocus} object or a list of such.}

\item{matchLoci}{A logical. If TRUE, the column names of \code{datamatrix} must be
found either within \code{names(loci)} or within the \code{name} entries of \code{loci}.
The column names of \code{datamatrix} are assumed to come in pairs with suffixes
".1" and ".2", e.g. "TH01.1", "TH01.2", etc. If FALSE (the default) it is
assumed that the \code{loci} correspond to the (pairs of) columns in
\code{datamatrix} sequentially.}
}
\value{
A \code{ped} object, or a list of such.
}
\description{
\code{Familias} is a widely used software for forensic kinship computations, which
also features an accompanying R package (also called \code{Familias}). The
function documented here converts pedigrees and marker data from the R
version of \code{Familias} to \code{\link[pedtools:ped]{pedtools::ped()}} format, used by \code{forrel}. This may
be of interest for specialized computations not implemented in \code{Familias},
e.g. conditional simulations. Note: For importing ".fam" files created by
(the stand-alone) \code{Familias}, see \code{\link[=readFam]{readFam()}}.
}
\details{
The \code{Familias} program represents pedigrees and marker data in a way that
differs from the \code{ped} format in several ways, mostly because of the latter's
stricter definition of a \emph{pedigree}. A \code{ped} object always represent a
connected pedigree, and each member must have either 0 or 2 parents. None of
this is required by \code{FamiliasPedigree} objects. The conversion function
\code{Familias2ped} takes care of all potential differences: It converts each
\code{Familias} pedigree into a list of connected \code{ped} objects, adding missing
parents where needed.
}
\examples{

famPed = structure(
  list(id = c('mother', 'daughter', 'AF'),
       findex = c(0, 3, 0),
       mindex = c(0, 1, 0),
       sex = c('female', 'female', 'male')),
  class = "FamiliasPedigree")

datamatrix = data.frame(
  M1.1 = c(NA, 8, NA),
  M1.2 = c(NA, 9.3, NA),
  row.names = famPed$id)

famLoc = structure(
  list(locusname = "M1",
       alleles = c("8" = 0.2, "9" = 0.5, "9.3" = 0.3)),
  class = "FamiliasLocus")

Familias2ped(famPed, datamatrix, loci = famLoc, matchLoci = TRUE)


}
\references{
Familias is freely available from \url{https://familias.name}.
}
\seealso{
\code{\link[=readFam]{readFam()}}
}
\author{
Magnus Dehli Vigeland, Thore Egeland
}
