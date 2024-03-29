% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boin.R
\name{Boin-class}
\alias{Boin-class}
\alias{Boin}
\title{The BOIN design as a Cumulative-Cohort Design (CCD) subclass}
\description{
The BOIN design as a Cumulative-Cohort Design (CCD) subclass

The BOIN design as a Cumulative-Cohort Design (CCD) subclass
}
\details{
TODO: Provide references, and note clearly the choice of defaults.
Also provide citations to BOIN paper(s).
}
\examples{

## ------------------------------------------------
## Method `Boin$new`
## ------------------------------------------------

# TODO
}
\section{Super classes}{
\code{\link[precautionary:Cpe]{precautionary::Cpe}} -> \code{\link[precautionary:Ccd]{precautionary::Ccd}} -> \code{Boin}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Boin$new()}}
\item \href{#method-clone}{\code{Boin$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="J">}\href{../../precautionary/html/Cpe.html#method-J}{\code{precautionary::Cpe$J()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="bU">}\href{../../precautionary/html/Cpe.html#method-bU}{\code{precautionary::Cpe$bU()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="max_dose">}\href{../../precautionary/html/Cpe.html#method-max_dose}{\code{precautionary::Cpe$max_dose()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="path_array">}\href{../../precautionary/html/Cpe.html#method-path_array}{\code{precautionary::Cpe$path_array()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="path_matrix">}\href{../../precautionary/html/Cpe.html#method-path_matrix}{\code{precautionary::Cpe$path_matrix()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="path_probs">}\href{../../precautionary/html/Cpe.html#method-path_probs}{\code{precautionary::Cpe$path_probs()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="path_rx">}\href{../../precautionary/html/Cpe.html#method-path_rx}{\code{precautionary::Cpe$path_rx()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="report">}\href{../../precautionary/html/Cpe.html#method-report}{\code{precautionary::Cpe$report()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Cpe" data-id="trace_paths">}\href{../../precautionary/html/Cpe.html#method-trace_paths}{\code{precautionary::Cpe$trace_paths()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="precautionary" data-topic="Ccd" data-id="applied">}\href{../../precautionary/html/Ccd.html#method-applied}{\code{precautionary::Ccd$applied()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Boin$new(target, cohort_max, enroll_max)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target}}{Target toxicity rate}

\item{\code{cohort_max}}{Upper bound on dose-wise enrollment}

\item{\code{enroll_max}}{Upper bound on total enrollment}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Create a new \code{Boin} object.
}

\subsection{Returns}{
A Boin object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# TODO
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Boin$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
