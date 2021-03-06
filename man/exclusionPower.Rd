% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exclusionPower.R
\name{exclusionPower}
\alias{exclusionPower}
\title{Power of exclusion}
\usage{
exclusionPower(
  claimPed,
  truePed,
  ids,
  markers = NULL,
  source = "claim",
  disableMutations = NA,
  exactMaxL = Inf,
  nsim = 1000,
  seed = NULL,
  alleles = NULL,
  afreq = NULL,
  knownGenotypes = NULL,
  Xchrom = FALSE,
  plot = FALSE,
  plotMarkers = NULL,
  verbose = TRUE
)
}
\arguments{
\item{claimPed}{A \code{ped} object (or a list of such), describing the claimed
relationship. If a list, the sets of ID labels must be disjoint, that is,
all ID labels must be unique.}

\item{truePed}{A \code{ped} object (or a list of such), describing the true
relationship. ID labels must be consistent with \code{claimPed}.}

\item{ids}{Individuals available for genotyping.}

\item{markers}{A vector indicating the names or indices of markers attached
to the source pedigree. If NULL (default), then all markers attached to the
source pedigree are used. If \code{alleles} or \code{afreq} is non-NULL, then this
parameter is ignored.}

\item{source}{Either "claim" (default) or "true", deciding which pedigree is
used as source for marker data.}

\item{disableMutations}{This parameter determines how mutation models are
treated. Possible values are as follows:
\itemize{
\item \code{NA} (the default): Mutations are disabled only for those markers whose
known genotypes are compatible with both \code{claimPed} and \code{truePed}. This is
determined by temporarily removing all mutation models and checking which
markers have nonzero likelihood in both alternatives.
\item \code{TRUE}: Mutations are disabled for all markers.
\item \code{FALSE}: No action is done to disable mutations.
\item A vector containing the names or indices of those markers for which
mutations should be disabled.
}}

\item{exactMaxL}{A positive integer, or \code{Inf} (default). Exact EPs are
calculated for markers whose number of alleles is less or equal to
\code{exactMaxL}; remaining markers are handled by simulation.}

\item{nsim}{A positive integer; the number of simulations used for markers
whose number of alleles exceeds \code{exactMaxL}.}

\item{seed}{An integer seed for the random number generator (optional).}

\item{alleles, afreq, Xchrom}{If these are given, they are used (together with
\code{knownGenotypes}) to create a marker object on the fly.}

\item{knownGenotypes}{A list of triplets \verb{(a, b, c)}, indicating that
individual \code{a} has genotype \code{b/c}. Ignored unless \code{alleles} or \code{afreq} is
non-NULL.}

\item{plot}{Either a logical or the character "plotOnly". If the latter, a
plot is drawn, but no further computations are done.}

\item{plotMarkers}{A vector of marker names or indices whose genotypes are to
be included in the plot.}

\item{verbose}{A logical.}
}
\value{
If \code{plot = "plotOnly"}, the function returns NULL after producing the
plot.

Otherwise, the function returns an \code{EPresult} object, which is essentially
a list with the following entries:
\itemize{
\item \code{EPperMarker}: A numeric vector containing the exclusion power of each
marker. If the known genotypes of a marker are incompatible with the true
pedigree, the corresponding entry is \code{NA}.
\item \code{EPtotal}: The total exclusion power, computed as \code{1 - prod(1 -   EPperMarker, na.rm = TRUE)}.
\item \code{expectedMismatch}: The expected number of markers giving exclusion,
computed as \code{sum(EPperMarker, na.rm = TRUE)}.
\item \code{distribMismatch}: The probability distribution of the number of markers
giving exclusion. This is given as a numeric vector of length \code{n+1}, where
\code{n} is the number of nonzero elements of \code{EPperMarker}. The vector has
names \code{0:n}.
\item \code{time}: The total computation time.
\item \code{params}: A list containing the (processed) parameters \code{ids}, \code{markers}
and \code{disableMutations}.
}
}
\description{
Computes the power (of a single marker, or for a collection of markers) of
excluding a claimed relationship, given the true relationship.
}
\details{
This function implements the formula for exclusion power as defined and
discussed in (Egeland et al., 2014).

It should be noted that \code{claimPed} and \code{truePed} may be any (lists of)
pedigrees, as long as they both contain the individuals specified by \code{ids}.
In particular, either alternative may have inbred founders (with the same or
different coefficients), but this must be set individually for each.
}
\examples{

############################################
### A standard case paternity case:
### Compute the power of exclusion when the claimed father is in fact
### unrelated to the child.
############################################

# Claim: Individual 1 is the father of individual 3
claim = nuclearPed(nch = 1, sex = 2)

# Truth: 1 and 3 are unrelated
true = list(singleton(id = 1), singleton(id = 3, sex = 2))

# Attach two markers
m1 = marker(claim, alleles = 1:2)
m2 = marker(claim, alleles = 1:3)
claim = setMarkers(claim, list(m1, m2))

# Compute EP when father and child is available for genotyping
exclusionPower(claim, true, ids = c(1,3))

# Suppose child is already genotyped
genotype(claim, marker = 1, id = 3) = c(1, 1)
genotype(claim, marker = 2, id = 3) = c(1, 1)

exclusionPower(claim, true, ids = 1)


############################################
### Two females claim to be mother and daughter, but are in reality sisters.
### We compute the power of various markers to reject the claim.
############################################

mother_daughter = nuclearPed(1, sex = 2)
sisters = relabel(nuclearPed(2, sex = c(2, 2)), c(101, 102, 2, 3))
ids = 2:3

# SNP with MAF = 0.1:
PE1 = exclusionPower(claimPed = mother_daughter, truePed = sisters,
                     ids = ids, alleles = 2, afreq = c(0.9, 0.1))

# Tetra-allelic marker with one major allele:
PE2 = exclusionPower(claimPed = mother_daughter, truePed = sisters,
                     ids = ids, alleles = 4, afreq = c(0.7, 0.1, 0.1, 0.1))

stopifnot(all.equal(c(PE1$EPtotal, PE2$EPtotal), c(0.00405, 0.03090)))

### How does the power change if the true pedigree is inbred?
sisters_LOOP = addParents(sisters, 101, father = 201, mother = 202)
sisters_LOOP = addParents(sisters_LOOP, 102, father = 201, mother = 203)


# SNP with MAF = 0.1:
PE3 = exclusionPower(claimPed = mother_daughter, truePed = sisters_LOOP,
                     ids = ids, alleles = 2, afreq = c(0.9, 0.1))

stopifnot(all.equal(PE3$EPtotal, 0.00765))

}
\references{
T. Egeland, N. Pinto and M.D. Vigeland, \emph{A general approach to
power calculation for relationship testing.} Forensic Science
International: Genetics 9 (2014): 186-190.
\doi{10.1016/j.fsigen.2013.05.001}
}
\author{
Magnus Dehli Vigeland
}
