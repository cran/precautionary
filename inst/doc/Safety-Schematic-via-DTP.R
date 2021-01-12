## ----setup, include=FALSE-----------------------------------------------------
old <- options(rmarkdown.html_vignette.check_title = FALSE) # suppress un-needed warning
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.height = 4, fig.width = 6)
knitr::opts_knit$set(eval.after = "fig.cap")

library(precautionary)
library(knitr)
library(kableExtra)
library(dplyr)
library(dtpcrm)
library(latticeExtra)

options(ordinalizer = NULL) # vignette presumes this is not already set

# Echo source focused on package functionality by redacting calls to kable()
# see https://bookdown.org/yihui/rmarkdown-cookbook/hook-hide.html:
local({
  hook_source <- knitr::knit_hooks$get('source')
  knitr::knit_hooks$set(source = function(x, options) {
    x <- gsub(" -> etc.", "", x)
    x <- x[!grepl("^etc. %>%", x)] # strip lines starting "etc. %>% ..."
    x <- x[!grepl("\\(etc.,", x)]  # strip lines with fun(etc., ...)
    hook_source(x, options)
  })
})

## ----clinical-params----------------------------------------------------------
number.doses <- 7
start.dose.level <- 3
max.sample.size <- 21
target.DLT <- 0.2
cohort.size <- 3

## ----model-spec-params--------------------------------------------------------
prior.DLT <- c(0.03, 0.07, 0.12, 0.20, 0.30, 0.40, 0.52)
prior.var <- 0.75

## ----early-stopping-----------------------------------------------------------
stop_func <- function(x) {
  y <- stop_for_excess_toxicity_empiric(x,
                                        tox_lim = target.DLT + 0.1,
                                        prob_cert = 0.72,
                                        dose = 1,
                                        nsamps = 100000)
  if(y$stop){
    x <- y
  } else {
    x <- stop_for_consensus_reached(x, req_at_mtd = 12)
  }
}

## ----compute-dtp, eval=FALSE--------------------------------------------------
#  t0 <- proc.time()
#  viola_dtp <- calculate_dtps(next_dose = start.dose.level,
#                              cohort_sizes = rep(cohort.size, 7),
#                              dose_func = applied_crm,
#                              prior = prior.DLT,
#                              target = target.DLT,
#                              stop_func = stop_func,
#                              scale = sqrt(prior.var),
#                              no_skip_esc = TRUE,
#                              no_skip_deesc = FALSE,
#                              global_coherent_esc = TRUE
#                              )
#  save(viola_dtp, file="../data/viola_dtp.rda")
#  proc.time() - t0 # ~17 minutes on a 2.6GHz i7

## ----load-cached-dtp, echo=FALSE----------------------------------------------
data(viola_dtp) # Read cached viola_dtp from disk

## ----degenerates--------------------------------------------------------------
knitr::kable(viola_dtp[1000:1010,])

## ----degeneracy---------------------------------------------------------------
viola_paths <- as.matrix(unique(viola_dtp))
nrow(viola_paths)

## ----T-constructed------------------------------------------------------------
I <- outer(viola_paths[,paste0("D",0:6)], 1:7, FUN = "==")
I[I] <- 1
I[!I] <- NA
# Now I[j,c,d] is an indicator array that we may use
# in the manner of a bitmask, to select the toxicities
# into their proper positions within T[j,c,d]:
T <- I * outer(viola_paths[,paste0("T",1:7)], rep(1,7))
dimnames(T)[[2]] <- paste0("C",1:7) # best labeled as Cohorts
dimnames(T)[[3]] <- paste0("X",1:7) # label as doses X1..X7
dim(T)

## ----b-via-T------------------------------------------------------------------
b <- apply(log(choose(3, T)), MARGIN = 1, FUN = sum, na.rm = TRUE)
length(b)

## ----Y-soeasy-----------------------------------------------------------------
Y <- apply(T, MARGIN = c(1,3), FUN = sum, na.rm = TRUE)
Z <- apply(3-T, MARGIN = c(1,3), FUN = sum, na.rm = TRUE)
U <- cbind(Y, Z)
dim(U)

## ----pi-per-skeleton----------------------------------------------------------
log_p <- log(prior.DLT)
log_q <- log(1 - prior.DLT)
log_pi <- b + U %*% c(log_p, log_q)
sum(exp(log_pi)) # check probabilities sum to 1

## ----VIOLA-dose-scaling-------------------------------------------------------
viola_dosing <-
  data.frame(Label = -2:4 # VIOLA dose designations -2,...,4
            ,Level =  1:7  # We use integer dose indexes here
            ,Dose_mg = c(0, 2.5, 5, 10, 15, 25, 35)
            ,skeleton = prior.DLT
             )
viola_dosing <- viola_dosing %>%
  mutate(log_dose = log(Dose_mg)
        ,sqrt_dose = sqrt(Dose_mg)
        )

plot(Level ~ sqrt_dose, data = viola_dosing, type = 'b'
     , ylab = "Dose Level"
     , xlab = expression(sqrt(Dose[mg]))
     , las = 1
     )

## ----skeleton-as-expectation, fig.cap="CRM skeleton probabilities vs $\\sqrt{\\mathrm{Dose}}$, with superimposed normal distribution having median $\\sqrt{33}$ and standard deviation 3."----
plot(skeleton ~ sqrt_dose, data = viola_dosing
     , type = 'b'
     , xlab = expression(sqrt(Dose[mg]))
     , ylim = c(0, 0.55)
     , las = 1)
x <- seq(0, 6, 0.1)
y <- pnorm(x, mean = sqrt(33), sd = 3)
lines(y ~ x, lty = 2)

## ----skeleton-log-match, fig.cap=paste0("CRM skeleton probabilities vs $\\log \\mathrm{Dose}$, with superimposed lognormal distribution having median $\\mu = ", exp(mu_opt), "$mg and $\\sigma = ", sigma_opt, "$.")----
plot(skeleton ~ log_dose, data = subset(viola_dosing, is.finite(log_dose))
     , type = 'b'
     , xlab = expression(log(Dose[mg]))
     , xlim = c(0.5, 3.6)
     , ylim = c(0, 0.55)
     , las = 1)

sse <- function(mu_sigma, dose = viola_dosing$Dose_mg){
  mu <- mu_sigma[1]
  sigma <- mu_sigma[2]
  lnorm.DLT <- plnorm(dose, meanlog = mu, sdlog = sigma)
  sse <- sum((lnorm.DLT - prior.DLT)[dose > 0]^2)
}
fit <- optim(par = c(mu=log(38), sigma=1.8), fn = sse)
mu_opt <- log(round(exp(fit$par['mu']), 1))
sigma_opt <- round(fit$par['sigma'], 1)
x <- seq(0.5, 3.55, 0.05)
y <- plnorm(exp(x), meanlog = mu_opt, sdlog = sigma_opt)
lines(y ~ x, lty = 2)

## ----fatalities-vs-kappa------------------------------------------------------
f <- function(kappa, mu=mu_opt, sigma=sigma_opt, X=viola_dosing$Dose_mg) {
  # For ease of use, this function is VECTORIZED over kappa,
  # generally returning a length(X) x length(kappa) MATRIX
  # which can validly right-multiply Y in t(pi) %*% Y %*% f.
  # (In the special case length(kappa)==1, a vector is returned.)
  gr5 <- plnorm(outer(X, exp(-2*kappa)), meanlog = mu, sdlog = sigma)
  dlt <- plnorm(outer(X, exp( 0*kappa)), meanlog = mu, sdlog = sigma)
  ff <- gr5/dlt
  ff[is.nan(ff)] <- 0 # replace NaN caused by X[1]==0
  dimnames(ff) <- list(
    dose = paste0("X", seq_along(X)),
    kappa = round(kappa, 3)
  )
  if (length(kappa)==1)
    return(as.vector(ff))
  ff
}

kappa <- seq(0.2, 1.2, 0.02)
F <- f(kappa = kappa)
Ef <- t(exp(log_pi)) %*% Y %*% F
Ef <- t(Ef) # transpose to a column-vector
plot(Ef ~ kappa, type = 'l'
     , xlab = expression(kappa)
     , ylab = "Expected Number of Fatal DLTs"
     , las=1)

## ----log_pi-------------------------------------------------------------------
log_pi <- function(mu, sigma){
  p <- plnorm(viola_dosing$Dose_mg, meanlog = mu, sdlog = sigma)
  log_pq <- c(log(p), log(1-p))
  log_pi = b + U %*% pmax(log_pq, -500) # clamping -Inf to -500 avoids NaN's  
}

## ----minimax------------------------------------------------------------------
mu_minimax <- log(viola_dosing$Dose_mg)[3] # dose 3 is 2nd-highest *nonzero* dose

## ----delta--------------------------------------------------------------------
delta <- mean(diff(log(viola_dosing$Dose_mg)[2:4]))

## ----kappa-delta-plane--------------------------------------------------------
focustab <- CJ(mu = mu_minimax
              , K = seq(0.5, 1.65 , 0.01) # K = kappa / sigma
              , sigma = delta/seq(0.4, 2.2, 0.01)
              )
focustab[, kappa := K * sigma]
focustab[, value := t(exp(log_pi(mu,sigma))) %*%
             Y %*% f(kappa=kappa, mu=mu, sigma=sigma)
       , by = .(mu,sigma)]

## ----contourplot, echo=FALSE, fig.asp=1, fig.cap="Ex ante expectation of fatalities in the VIOLA trial, under the analytical set-up of Figure 3 in @norris_what_2020. Therapeutic index $\\kappa/\\sigma$ gauges the target drug's aptness for safe-and-effective 1-size-fits-all dosing in the trial's combination regimen. Signal-to-noise index $\\delta/\\sigma$ governs the informativeness of the dose-escalation process. This analysis assumes a log-normally distributed $\\MTDi$, with median parameter set (on minimax grounds) at VIOLA dose level 3. Plotted values of $\\delta/\\sigma$ are based on $\\delta=\\log(2)$, selected according to the $2\\times$ VIOLA dose-level spacing around this dose level."----
contourplot(value ~ K + (delta/sigma)
          , data = focustab
          , at = seq(0.0, 2.0, 0.1)
          , par.strip.text = list(cex=0.7)
          , xlab = list(expression(kappa / sigma), rot=0)
          , ylab = list(expression(frac(delta, sigma)), rot=0)
          , scale = list(cex=0.7, x=list(at=c(0.5, 1, 1.5)))
          , label.style = "align"
          , labels = list(cex=0.5)
          , aspect = 1
          , col.regions = hcl.colors(10)
            )

## ----echo=FALSE, results='hide'-----------------------------------------------
options(old) # restore user's original options before finishing, per CRAN

## ----bib, include=FALSE, cache=FALSE------------------------------------------
# Create a bib file for packages cited in this paper
knitr::write_bib(c('dtpcrm'), file = 'packages.bib')

