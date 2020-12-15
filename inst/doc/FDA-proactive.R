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

## ----design-------------------------------------------------------------------
library(precautionary)
design <- get_three_plus_three(num_doses = 6)
options(dose_levels = c(2, 6, 20, 60, 180, 400)) # ng/kg/week

## ----hyperprior, fig.cap="Multiple samples from a *hyperprior* over the distribution of MTD$_i$."----
mtdi_gen <- hyper_mtdi_lognormal(CV = 0.5
                                ,median_mtd = 180
                                ,median_sdlog = 0.6
                                ,units="ng/kg/week"
                                )
plot(mtdi_gen, n=100, col=adjustcolor("red", alpha=0.25))

## ----ordinalizer--------------------------------------------------------------
options(ordinalizer = function(MTDi, r0 = 1.5) {
  MTDi * r0 ^ c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
})

## ----simulate-----------------------------------------------------------------
set.seed(2020)
design %>% simulate_trials(
  num_sims = 100
, true_prob_tox = mtdi_gen
) %>% extend(target_mcse = 0.1) -> SIMS

## ----summarize----------------------------------------------------------------
summary(SIMS,r0 = 2)$safety %>%
  safety_kable() 

## ----tabulate-safety----------------------------------------------------------
r0 <- seq(1.2, 2.0, 0.2); names(r0) <- format(r0, 2)
safetytabs <- sapply(r0, FUN = function(.) summary(SIMS, r0=.)$safety
                     , simplify = "array", USE.NAMES = TRUE)
cbind(data.table(`$r_0$` = r0), t(safetytabs[1,,])) %>%
  kable(digits=2) %>% add_header_above(
    c(" "=1, "Expected counts by toxicity grade"=6, " "=1))

## ----focus-fatality-----------------------------------------------------------
cbind(data.table(`$r_0$` = r0), t(safetytabs[,'Gr5',])) %>%
  kable(digits=2,
        col.names = c('$r_0$',
                      'Expected fatal toxicities',
                      'MCSE')) %>%
  kable_styling(full_width = FALSE, position = "left")

## ----echo=FALSE, results='hide'-----------------------------------------------
options(old) # restore user's original options before finishing, per CRAN

