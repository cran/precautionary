## ----setup, include=FALSE-----------------------------------------------------
old <- options(rmarkdown.html_vignette.check_title = FALSE) # suppress un-needed warning
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.height = 4, fig.width = 6)

library(precautionary)
library(knitr)
library(kableExtra)

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

## ----basic-sim----------------------------------------------------------------
library(escalation)
# Design a 3 + 3 trial, with 5 pre-specified doses to escalate through
design <- get_three_plus_three(num_doses = 5)
# Posit a scenario where these 5 doses cause dose-limiting toxicities
# (DLTs) in 12%, 27%, etc. of the population:
scenario <- c(0.12, 0.27, 0.44, 0.53, 0.57)
design %>% simulate_trials( # Feed the design to the simulator ...
  num_sims = 100            # run 100 simulated trials
, true_prob_tox = scenario  # under the chosen scenario,
) -> sims                   # and store the results.

summary(sims) -> etc. # summarize simulation results
etc. %>%
  kable(caption="Simulation summary as produced by package `escalation`.")

## -----------------------------------------------------------------------------
library(precautionary)
mtdi_dist <- mtdi_lognormal(CV = 2          # coefficient of variation
                           ,median = 5      # median DLT threshold
                           ,units = "mg/kg" # real doses have units!
                           )

## -----------------------------------------------------------------------------
options(dose_levels = c(0.5, 1, 2, 4, 6)) # specify actual dosing

## -----------------------------------------------------------------------------
plot(mtdi_dist)

## -----------------------------------------------------------------------------
probs <- mtdi_dist@dist$cdf(getOption('dose_levels'))
names(probs) <- paste(getOption('dose_levels'), mtdi_dist@units)
t(probs) %>% kable(digits = 4)

## ----MTDi-regime--------------------------------------------------------------
design %>% simulate_trials(
  num_sims = 100
, true_prob_tox = mtdi_dist # pull tox probs from a MODEL, not thin air
) -> SIMS

summary(SIMS) -> etc.
etc. %>% kable(caption = "A simulation summary under the MTD$_i$ regime of package `precautionary`.")

## ---- echo=FALSE--------------------------------------------------------------


## -----------------------------------------------------------------------------
tox_threshold_scaling <- function(MTDi, r0) {
  MTDi * r0 ^ c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
}
summary(SIMS
       ,ordinalizer = tox_threshold_scaling
       ,r0 = 2      # supply a value for ordinalizer's r0 parameter
       )$safety %>% # select the 'safety' component of the summary
  safety_kable()

## ---- eval=FALSE--------------------------------------------------------------
#  tox_threshold_scaling <-
#    function(MTDi   # An ordinalizer is a function of a dose threshold,
#            ,r0 = 2 # and in general has additional parameters as well.
#            ) {
#      # An ordinalizer assumes we start with a binary toxicity notion,
#      # and maps that to a *graded* notion of toxicity by means of a
#      # transformation in 'dose-space'.
#      # Assuming the default value r0 = 2 provided in its definition,
#      # this ordinalizer says that an individual whose dose threshold
#      # for the binary toxicity is MTDi has thresholds ...
#      c(Gr1 = MTDi / r0^2 # at MTDi/4 for Gr1,
#       ,Gr2 = MTDi / r0   # at MTDi/2 for Gr2,
#       ,Gr3 = MTDi        # at MTDi for Gr3 ('tox' is defined as Gr3+),
#       ,Gr4 = MTDi * r0   # at 2*MTDi for Gr4,
#       ,Gr5 = MTDi * r0^2 # at 4*MTDi for Gr5.
#      )
#    }

## ---- eval=FALSE--------------------------------------------------------------
#  function(MTDi3, r12, r4, r5) {
#    c(Gr1 = MTDi3 / r12^2
#     ,Gr2 = MTDi3 / r12
#     ,Gr3 = MTDi3
#     ,Gr4 = MTDi3 * r4
#     ,Gr5 = MTDi3 * r4*r5
#     )
#  }

## ----Raleigh-distribution, echo=FALSE, message=FALSE--------------------------
# Plot the Raleigh density for a suitable range of sigma_CV values
sigmas <- c(0.1, 0.25, 0.5, 1)
curves <- as.data.table(expand.grid(x = seq(0, 2, 0.01)
                                   ,mode = sigmas
                                   )
                        )
curves[, pdf := 2*(x/mode^2)*dchisq((x/mode)^2, df=2), by = mode]
plot(pdf ~ x, data = curves[mode == sigmas[1]], type = 'l'
     , main = "Probability density function of the Rayleigh distribution")
for(sigma in sigmas[-1])
  lines(pdf ~ x, data = curves[mode == sigma])
text(x=0.3, y=5, expression(sigma == 0.1))
text(x=0.5, y=2.25, expression(sigma == 0.25))
text(x=0.75, y=1.4, expression(sigma == 0.5))
text(x=1.5, y=0.9, expression(sigma == 1))

## ----hyperprior, fig.cap="Multiple samples from a *hyperprior* over the distribution of MTD$_i$. Consider the implications for the customary practice of pulling 'true toxicity probabilities' out of thin air!"----
mtdi_gen <- hyper_mtdi_lognormal(CV = 1
                                ,median_mtd = 5
                                ,median_sdlog = 0.5 # this is new
                                ,units="mg/kg"
                                )
plot(mtdi_gen, n=100, col=adjustcolor("red", alpha=0.25))

## ----dang---------------------------------------------------------------------
design %>% simulate_trials(
  num_sims = 400
, true_prob_tox = mtdi_gen # pull tox probs from MANY models
) -> HYPERSIMS

# As a convenience, package 'precautionary' lets you set the
# ordinalizer as an *option* so that you don't have to keep
# specifying it as an argument to summary(). By providing a
# default setting for the r0 parameter in the ordinalizer
# definition, we avoid having to keep specifying that, too.
options(ordinalizer = function(MTDi, r0 = 1.5) {
  MTDi * r0 ^ c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
})

summary(HYPERSIMS)$safety -> etc.
etc. %>% safety_kable()

## -----------------------------------------------------------------------------
r0 <- c(1.25, 1.5, 1.75, 2.0)
rbind(summary(HYPERSIMS, r0=r0[1])$safety[1,]
     ,summary(HYPERSIMS, r0=r0[2])$safety[1,]
     ,summary(HYPERSIMS, r0=r0[3])$safety[1,]
     ,summary(HYPERSIMS, r0=r0[4])$safety[1,]
     ) -> safety
cbind(data.table(`$r_0$` = r0), safety) %>% kable(digits=2) %>%
  add_header_above(c(" "=1, "Expected counts by toxicity grade"=6, " "=1))

## ----hyper-CRM----------------------------------------------------------------
crm_design <- get_dfcrm(skeleton = scenario, target = 0.25) %>%
  stop_at_n(n = 24)
crm_design %>% simulate_trials(
  num_sims = 200
, true_prob_tox = mtdi_gen
) -> CRM_HYPERSIMS

summary(CRM_HYPERSIMS)$safety -> etc.
etc. %>% safety_kable()

## ----hyper-BOIN---------------------------------------------------------------
boin_design <- get_boin(num_doses = 5, target = 0.25) %>%
  stop_at_n(n = 24)
boin_design %>% simulate_trials(
  num_sims = 200
, true_prob_tox = mtdi_gen
) -> BOIN_HYPERSIMS

summary(BOIN_HYPERSIMS)$safety -> etc.
etc. %>% safety_kable()

## ----echo=FALSE, results='hide'-----------------------------------------------
options(old) # restore user's original options before finishing, per CRAN

