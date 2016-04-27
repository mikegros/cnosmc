\name{cno_smc}
\alias{cno_smc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to perform SMC on the CNO protein network}
\description{TO DO: Add more of a full description later}
\usage{
cno_smc(n_samples, data, model, init_links = 1, p_link = 0.9, n_mh = 5, split_inhib = FALSE, n_cores = 1, diagnostics = F)
}
\arguments{
  \item{n_samples}{scalar number of SMC samples generated}
  \item{data}{data object from CNOFuzzy network. Contains simulated data, input settings, etc.}
  \item{model}{model object from CNOFuzzy network. Contains main values for showing the network structure}
  \item{init_links}{scalar or vector of indices of active links at the start of the sampling. Must contain one link to a stimulus node and one to an observed node and be a completed path.}
  \item{p_link}{scalar Prior probability that a link is active}
  \item{n_mh}{scalar number of Metropolis-Hastings jittering steps for preventing particle degeneracy in the SMC.}
  \item{split_inhib}{logical indicator of whether Gstring parameters should be split for each setting combination for the inhibitors or not. Allows model structure to be identified conditional on the input settings of the inhibitors. This is valuable in dealing with identifiability of links.}
  \item{n_cores}{scalar number of parallel cores to use for sampling through the SMC particles. n_cores = 1 runs in serial. NOTE (MJG): PARALLEL HAS NOT YET BEEN TESTED!}
  \item{diagnostics}{logical indicator of whether diagnostic information should be saved or printed. NOTE(MJG): Not yet tested and will give a lot of junk if TRUE for now}
}
\details{TO DO: Add details}
\value{list of matrices containing the SMC samples for the transfer function parameters on the links as well as for the indicator of whether the link is active.}
\references{TO DO: Add references}
\author{Mike Grosskopf}
\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Sequential Monte Carlo}
