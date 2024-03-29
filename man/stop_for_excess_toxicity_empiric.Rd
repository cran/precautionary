% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dtp.R
\name{stop_for_excess_toxicity_empiric}
\alias{stop_for_excess_toxicity_empiric}
\title{A supremely faster version of a function from \code{dtpcrm} v0.1.1}
\usage{
stop_for_excess_toxicity_empiric(
  x,
  tox_lim,
  prob_cert,
  dose = 1,
  suppress_dose = TRUE
)
}
\arguments{
\item{x}{A object of class \code{mtd}}

\item{tox_lim}{Scalar upper threshold on estimated toxicity rate}

\item{prob_cert}{Confidence level for threshold exceedance}

\item{dose}{Integer scalar, the dose being considered}

\item{suppress_dose}{Logical; if TRUE the MTD is set to \code{NA} when
trial stop is recommended.}
}
\value{
The \code{mtd} object x, with stop decision annotated
}
\description{
Originally, the sampling in \code{stats::rnorm()} (see inline comments in code)
consumed 53\% of run-time in a 6-cohort VIOLA DTP. After this change, it
doesn't even show up! More importantly, the consumption is now dominated
by (at 75\%) by the objective function 'f' in integrate().
}
\author{
Adapted by David C. Norris from original \code{dtpcrm::stop_for_excess_toxicity_empiric}
}
