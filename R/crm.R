## Introducing an R6 class for CRM models, capable of supporting
## efficient caching maneuvers.

## TODO: Organize a class hierarchy of model types (empiric, logistic),
##       exploiting inheritance to clarify what is special vs shared.

#' @name Crm-class
#' @title An R6 class encapsulating CRM models
#'
#' This class wraps the functionality of package \CRANpkg{dfcrm}, adding
#' efficient Rust implementations of some numerical routines.
#'
#' @details
#' Syntactically, the method chaining supported by R6 classes makes the
#' invocation of CRM models more transparent. The mutability conferred
#' by reference semantics enables memoization (caching) of results, which
#' can speed up DTP calculations significantly.
#'
#' Presently, this class supports only the 'empiric' (aka 'power') model.
#' But it is hoped that inheritance will assist in rendering other models
#' implemented in package \CRANpkg{dfcrm} clearly, with code reuse.
#' @importFrom R6 R6Class
#' @export
Crm <- R6Class("Crm",
               inherit = Cpe,
               public = list(
                 #' @details
                 #' Create a new `Crm` object.
                 #'
                 #' @param skeleton CRM skeleton
                 #' @param scale Sigma parameter of prior on beta parameter
                 #' @param target Target toxicity rate
                 #' @return A `Crm` object.
                 #'
                 #' @examples
                 #' # An example verbatim from dfcrm::crm()
                 #' prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
                 #' target <- 0.2
                 #' level <- c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2)
                 #' y     <- c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0)
                 #' s <- sqrt(1.34)
                 #' old <- dfcrm::crm(prior, target, y, level)
                 #' new <- Crm$new(skeleton = prior, target = target)$
                 #'          dontcache()$
                 #'          observe(level, y)$
                 #'          est(impl="rusti", abbrev=FALSE)
                 initialize = function(skeleton, scale = sqrt(1.34), target) {
                   private$ln_skel <- log(skeleton)
                   private$scale <- scale
                   private$target <- target
                   private$cache <- new.env(hash = TRUE, size = 10000L)
                   private$x <- private$o <- integer(length(private$ln_skel))
                 },
                 #' @details
                 #' Return number of prespecified doses
                 #' @note This specializes the superclass set/get method, consistent
                 #' with the non-mutable number of doses of CRM with fixed skeleton.
                 #'
                 #' @param D Included to match signature of superclass method.
                 #' It is an error to call this method with non-missing `D` parameter.
                 #' @return Length of CRM skeleton.
                 max_dose = function(D) {
                   if (missing(D))
                     return(length(private$ln_skel))
                   stop("Class 'Crm' has fixed number of doses determined by the skeleton.")
                 },
                 #' @details
                 #' Set or query CRM skeleton
                 #'
                 #' @param skeleton A numeric vector to set as the model skeleton.
                 #' @return Self (invisibly), unless `skeleton` is missing,
                 #' in which case the skeleton, a numeric vector, is returned.
                 skeleton = function(skeleton) {
                   if (missing(skeleton))
                     return(exp(private$ln_skel))
                   ## When *setting* the skeleton, we invalidate the cache:
                   private$cache <- new.env(hash = TRUE, size = 10000L)
                   private$skips <- 0
                   private$ln_skel <- log(skeleton)
                   super$max_dose(length(skeleton)) # TODO: Delete super's private$.max_dose?
                   invisible(self)
                 },
                 #' @details
                 #' Set private cache to NULL; useful for performance testing
                 #'
                 #' @return Self, invisibly
                 dontcache = function() {
                   private$cache <- NULL
                   private$skips <- NA
                   invisible(self)
                 },
                 #' @details
                 #' Report lifetime duty & performance statistics
                 #'
                 #' @param ... Optional arguments specifying columns to add to report
                 #' @return A named vector summarizing lifetime duty and performance
                 report = function(...) {
                   data.table(
                     pid = Sys.getpid(),
                     ...,
                     evals = private$evals,
                     skips = private$skips,
                     calc.ms = as.integer(1000*sum(private$user)),
                     '\u03bcs/calc' = as.integer(1000000*sum(private$user)/private$evals,3) # μs
                     )
                 },
                 #' @details
                 #' Set the stopping function
                 #'
                 #' @param sfunc A function taking `mtd` objects to `mtd` objects,
                 #' attaching suitable stopping information
                 #' @return Self, invisibly
                 stop_func = function(sfunc){
                   private$.stop_func <- sfunc
                   ## NB: Changing stop_func does NOT invalidate cache, since the
                   ##     dose recommmendation from $esc() gets cached *before*
                   ##     application of stopping rules occurs in $applied().
                   invisible(self)
                 },
                 #' @details
                 #' Set the `no_skip_esc` behavior
                 #'
                 #' @param tf An atomic logical value, TRUE or FALSE
                 #' @return Self, invisibly
                 no_skip_esc = function(tf){
                   stopifnot(is.logical(tf))
                   stopifnot(length(tf) == 1)
                   private$.no_skip_esc = tf
                   invisible(self)
                 },
                 #' @details
                 #' Set the `no_skip_deesc` behavior
                 #'
                 #' @param tf An atomic logical value, TRUE or FALSE
                 #' @return Self, invisibly
                 no_skip_deesc = function(tf){
                   stopifnot(is.logical(tf))
                   stopifnot(length(tf) == 1)
                   private$.no_skip_deesc = tf
                   invisible(self)
                 },
                 #' @details
                 #' Set the `global_coherent_esc` behavior
                 #'
                 #' @param tf An atomic logical value, TRUE or FALSE
                 #' @return Self, invisibly
                 global_coherent_esc = function(tf){
                   stopifnot(is.logical(tf))
                   stopifnot(length(tf) == 1)
                   private$.global_coherent_esc = tf
                   invisible(self)
                 },
                 #' @details
                 #' Set the required confidence level for escalation decisions
                 #'
                 #' @param conf A numeric confidence less than 1.0
                 #' @return Self, invisibly
                 conf_level = function(conf){
                   private$conf.level <- conf
                   invisible(self)
                 },
                 #' @details
                 #' Set patient-wise toxicity observations
                 #'
                 #' @param level A patient-wise vector of dose assignments
                 #' @param tox A patient-wise vector of 0/1 toxicity assessments
                 #' @return Self, invisibly
                 observe = function(level, tox){ # TODO: Preserve order for dfcrm's sake?
                   stopifnot(length(level) == length(tox))
                   stopifnot(all(tox %in% c(0,1)))
                   D <- self$max_dose()
                   stopifnot(all(level %in% 1:D))
                   level <- factor(level, levels=1:D)
                   self$tally(x = xtabs(tox ~ level),
                              o = xtabs(!tox ~ level))
                 },
                 #' @details
                 #' Set dose-wise toxicity observations
                 #'
                 #' @param x A dose-wise vector of toxicity counts
                 #' @param o A dose-wise vector of non-toxicity counts
                 #' @return Self, invisibly
                 tally = function(x, o){
                   D <- self$max_dose()
                   ## We expect 'x' to be a dosewise vector of exchangeable toxicity counts
                   stopifnot(length(x) == D)
                   stopifnot(is.numeric(x))
                   stopifnot(all(round(x) == x))
                   private$x <- as.integer(x)

                   ## And we expect a dosewise 'o'
                   stopifnot(length(o) == D)
                   ## that is either a numeric vector of exchangeable nontox counts ...
                   if (is.numeric(o)) {
                     private$o <- as.integer(o)
                     private$level <- c(rep(1:D, o),
                                        rep(1:D, x))
                     private$w <- c(rep(1, sum(o)),
                                    rep(0, sum(x)))
                   } else { # .. or a _list_ of weight vectors breaking exchangeability.
                     ## TODO: Consider a more functional construction of $w and $level:
                     private$w <- numeric(length(o) +
                                          sum(x)) # initialized to 0, thus encoding y
                     ix <- 0
                     for (level in 1:length(o)) {
                       len <- length(o[[level]])
                       private$w[ix + 1:len] <- o[[level]]
                       private$level[ix + 1:len] <- level
                       ix <- ix + len
                     }
                     private$o <- lapply(o, length)
                   }
                   invisible(self)
                 },
                 #' @details
                 #' Estimate the model
                 #'
                 #' @param impl A string choosing the low-level implementation to use.
                 #' Possible values include `"dfcrm"`, `"rusti"` and `"ruste"`.
                 #' @param abbrev Logical; if TRUE (the default), an abbreviated `mtd`
                 #' object is returned to save execution time. If FALSE, a complete object is
                 #' returned, suitable for regression testing against package \CRANpkg{dfcrm}.
                 #' @return An object of class `mtd` as per package \CRANpkg{dfcrm}
                 est = function(impl='rusti', abbrev=TRUE){
                   private$evals <- private$evals + 1
                   t0 <- proc.time()
                   scale <- private$scale
                   model <- "empiric"
                   method <- "bayes"
                   include <- seq_along(private$w)
                   pid <- include
                   switch(impl
                         ,dfcrm = {
                            ans <-
                              dfcrm::crm(prior = exp(private$ln_skel), target = private$target,
                                         tox=(private$w==0), level=private$level,
                                         n=length(private$w), dosename=NULL,
                                         include=include, pid=pid, conf.level=private$conf.level,
                                         method=method, model=model, intcpt=3,
                                         scale=scale, model.detail=TRUE, patient.detail=TRUE,
                                         var.est=TRUE)
                            ###private$user['dfcrm'] <- private$user['dfcrm'] +
                            ###  sum((proc.time() - t0)[c('user.self','user.child')])
                            return(ans)
                          }
                         ,rusti = {
                           ln_x1p <- private$ln_skel[private$level]
                           w <- private$w
                           m0 <- integrate(crmh  ,-Inf,Inf, ln_x1p, w, scale, abs.tol=0)[[1]]
                           m1 <- integrate(crmht ,-Inf,Inf, ln_x1p, w, scale, abs.tol=0)[[1]]
                           m2 <- integrate(crmht2,-Inf,Inf, ln_x1p, w, scale, abs.tol=0)[[1]]
                         }
                        ,ruste = {
                          ln_p <- private$ln_skel
                          tox <- private$x
                          nos <- private$o
                          m0 <- integrate(crmh_xo,-Inf,Inf,ln_p,tox,nos,scale,0,abs.tol=0)[[1]]
                          m1 <- integrate(crmh_xo,-Inf,Inf,ln_p,tox,nos,scale,1,abs.tol=0)[[1]]
                          m2 <- integrate(crmh_xo,-Inf,Inf,ln_p,tox,nos,scale,2,abs.tol=0)[[1]]
                        }
                       ,stop("must specify impl in Crm$est()")
                        )
                   private$user[impl] <- private$user[impl] +
                     sum((proc.time() - t0)[c('user.self','user.child')])
                   estimate <- m1/m0
                   post.var <- m2/m0 - estimate^2
                   prior <- exp(private$ln_skel)
                   ptox <- prior^exp(estimate)
                   ans <- list(prior = prior,       # used by stop_for_excess_toxicity_empiric
                               estimate = estimate, # used by stop_for_excess_toxicity_empiric
                               post.var = post.var, # used by stop_for_excess_toxicity_empiric
                               level = private$level, # used by stop_for_consensus_reached
                               mtd = which.min(abs(ptox - private$target)))
                   if (abbrev)
                     return(ans)
                   ## Otherwise we build out the whole return value of class "mtd"
                   crit <- qnorm(0.5 + private$conf.level/2)
                   lb <- estimate - crit*sqrt(post.var)
                   ub <- estimate + crit*sqrt(post.var)
                   ## TODO: Might 'within(.)' itself be costly? Try flattening
                   ##       this code to a series of straight component assignments.
                   ans <- within(ans, {
                     target <- private$target
                     conf.level <- private$conf.level
                     ptox <- ptox
                     ptoxL <- prior^exp(ub)
                     ptoxU <- prior^exp(lb)
                     tox <- 1.0*(private$w==0) # NB: dfcrm's "mtd" holds tox as double
                     dosename <- NULL
                     subset <- pid[include]
                     model <- model
                     prior.var <- scale^2
                     method <- method
                     include <- include
                     pid <- pid
                     model.detail <- TRUE
                     intcpt <- 3 # TODO: Un-hard-code this
                     ptoxL <- ptoxL
                     ptoxU <- ptoxU
                     patient.detail <- TRUE
                     tite <- FALSE
                     dosescaled <- prior # TODO: this is correct for EMPIRIC MODEL ONLY
                     var.est <- TRUE
                   })

                   ## Order components to enable testthat::expect_equal()
                   ans <- ans[c("prior","target","tox","level","dosename","subset",
                                "estimate","model","prior.var","post.var","method",
                                "mtd","include","pid","model.detail","intcpt",
                                "ptox","ptoxL","ptoxU","conf.level","patient.detail",
                                "tite","dosescaled","var.est")]
                   class(ans) <- "mtd"
                   return(ans)
                 }, #</est()>
                 #' @details
                 #' Return dose recommendation for given tox/no-tox tallies.
                 #'
                 #' This function caches results, which greatly saves computation time
                 #' in CPE -- yielding e.g. a 5x speedup for the VIOLA trial example.
                 #' @param x A dose-wise vector of toxicity counts
                 #' @param o A dose-wise vector of non-toxicity counts
                 #' @param last_dose The most recently given dose, as required to implement
                 #' the `global_coherent_esc=TRUE` behavior
                 #' @param max_dose Unused; included for compatibility with superclass method
                 #' of required `impl` parameter and optional `abbrev` flag.
                 #' @return An object of class `mtd` as per package \CRANpkg{dfcrm},
                 #' or possibly an abbreviated version of such object as returned by
                 #' method `Crm$est()`.
                 applied = function(x, o, last_dose = NA, max_dose = NULL){
                   if (!is.null(private$cache)) {
                     key <- paste(x, (x+o), sep='/', collapse='-') # human-readable to aid analysis
                     if (private$.global_coherent_esc)
                       key <- paste(key, last_dose, sep='@') # last_dose becomes relevant to lookup
                     if (!is.null(est <- get0(key, envir = private$cache))) {
                       private$skips <- private$skips + 1
                       return(est)
                     }
                   }
                   level <- which((x+o) > 0) # vector of levels that have been tried
                   est <- self$tally(x, o)$est()
                   if (private$.no_skip_esc) {
                     est$mtd <- min(est$mtd, max(level) + 1)
                   }
                   if (private$.no_skip_deesc) {
                     est$mtd <- max(est$mtd, min(level) - 1)
                   }
                   if (private$.global_coherent_esc) {
                     if(is.na(last_dose)) stop("last_dose required if global_coherent_esc = TRUE")
                     if (x[last_dose] > private$target * (x+o)[last_dose]) {
                       est$mtd <- min(est$mtd, last_dose)
                     }
                   }
                   if (!is.null(private$.stop_func)) {
                     est <- private$.stop_func(est)
                   }
                   ## Whereas 'dtpcrm' stopping functions recommend dose=NA when stopping
                   ## for excess toxicity, I find that a dose=0 rec in this case accords
                   ## nicely with the idiom of signalling stops with negative dose recs.
                   if(!is.na(est$stop) && est$stop && is.na(est$mtd))
                     est$mtd <- 0L
                   if (!is.null(private$cache))
                     assign(key, est, envir = private$cache)
                   return(est)
                 } #</applied>
               ), # </public>
               private = list(
                 ln_skel = NA
               , scale = NA
               , target = NA
               , .stop_func = NULL
               , .no_skip_esc = TRUE     # These would seem to be the
               , .no_skip_deesc = FALSE  # 'safest' defaults.
               , .global_coherent_esc = TRUE  # This sounds like a good thing.
               , x = NA
               , o = NA
               , w = NA # in general, this will encode y as well
               , level = NA # patient-wise dose level, such that ln_skel[l] aligns with w
               , obs = NaN
               , cache = NULL
               , evals = 0
               , skips = 0
               , user = c(dfcrm = 0.0,
                          rusti = 0.0,
                          ruste = 0.0)
               , conf.level = 0.90
               , path_list = list()
                 ## See Norris2020c for definitions of T, b, U, Y
               , T = NA # J*C*D array
               , b = NA # J-vector
               , U = NA # J*2D matrix
               , Y = NA # left half of U
               )
               )

## Parameters of objective functions are as follows:
## a: The single parameter of 1-parameter CRM. This function is
##    vectorized over this parameter, as required by 'integrate'.
## x: Vector of prior probabilities of tox at administered doses
## y: Vector of 0/1 indicators of toxicity, patient-wise
## w: A vector of weights?
## s: A (scalar) scale factor


#' A package-local (as-yet, unexported) test harness adapted from `dfcrm::crm()`.
#'
#' for various performance tuning experiments. The `impl` arg allows selection
#' of various alternative implementations:
#' - `'rusti'` substitutes integrands `crmh`, `crmht`, `crmht2` written in Rust
#' - `'dfcrm'` is the original as implemented in package `dfcrm`.
#' @param prior The CRM skeleton: dose-wise prior probabilities of toxicity
#' @param target Target toxicity rate
#' @param tox A patient-wise vector of toxicity counts
#' @param level A patient-wise vector of dose level assignments
#' @param n The number of patients enrolled
#' @param dosename Optional designators for the doses
#' @param include Index of patients to include
#' @param pid Vector of patient ID labels
#' @param conf.level Used to assign upper and lower bounds on predicted probability
#' of toxicity, which in turn may be referenced in escalation, deescalation and
#' stopping decisions.
#' @param method Estimation method:
#' @param model Presently, only the \sQuote{empiric} (or \sQuote{power}) model
#' has a Rust likelihood implementation.
#' @param intcpt Intercept for \sQuote{logistic} model
#' @param scale  Sigma parameter of prior on beta parameter
#' @param model.detail If FALSE, the model content of an `mtd` object will not
#' be displayed.  Default is TRUE.
#' @param patient.detail If FALSE, patient summary of an `mtd` object will not
#' be displayed.  Default is TRUE.
#' @param var.est If TRUE, variance of the estimate of the model parameter and
#' probability/confidence interval for the dose-toxicity curve will be computed
#' @param impl Switch between `'rusti'` and `'dfcrm'` implementations.
#' Currently the `'rusti'` option is implemented only for the Bayes method
#' of the empirical (\sQuote{power}) model. An experimental `'ruste'`
#' implementation is in the works.
#' @importFrom stats integrate optimize qnorm
#' @importFrom dfcrm crmhlgt crmhtlgt crmht2lgt
#' @importFrom dfcrm lcrm lcrmlgt
#' @author Adapted by David C. Norris, from Ken Cheung's \CRANpkg{dfcrm}
crm <- function(prior, target, tox, level, n=length(level),
                dosename=NULL, include=1:n, pid=1:n, conf.level=0.90,
                method="bayes", model="empiric", intcpt=3,
                scale=sqrt(1.34), model.detail=TRUE, patient.detail=TRUE, var.est=TRUE,
                impl=c("rusti","ruste","dfcrm")) { # implementation switch
  if (impl[1]=="dfcrm")
    return(dfcrm::crm(prior, target, tox, level, n, dosename, include, pid, conf.level,
                      method, model, intcpt, scale, model.detail, patient.detail, var.est))
  if (method != "bayes")
    warning("CRM method '", method, "' has no '", impl, "' implementation as yet.")
  if (model != "empiric")
    warning("The '", model, "' model has no '", impl, "' implementation as yet.")
  ## NB: If we reach this point without warnings, then we're estimating
  ##     the empiric model using Bayes method. This is the test case for
  ##     my alternative implementations.
  y1p <- tox[include]
  y1p <- as.integer(y1p) # Rust methods require this
  w1p <- rep(1,length(include))
  w1p[y1p == 1] <- 0; # encode y in w1p
  if (model=="empiric") {
    dosescaled <- prior

    x1p <- prior[level[include]]
    ln_x1p <- log(x1p) # Rust integrands {crmh,crmht,crmht2} take a log-skeleton to spare txops
    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrm,c(-10,10),x1p,y1p,w1p,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2,-Inf,Inf,ln_x1p,w1p,500,abs.tol=0)[[1]] / integrate(crmh,-Inf,Inf,x1p,y1p,w1p,500,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      switch(impl[1],
             rusti = {
               den <- integrate(crmh,-Inf,Inf,ln_x1p,w1p,scale,abs.tol=0)[[1]]
               est <- integrate(crmht,-Inf,Inf,ln_x1p,w1p,scale,abs.tol=0)[[1]] / den
               if (var.est)
                 e2 <- integrate(crmht2,-Inf,Inf,ln_x1p,w1p,scale,abs.tol=0)[[1]] / den
             },
             ruste = stop("impl='ruste' is inefficient unless invoked from class Crm"),
             stop(paste("impl =", impl, "not recognized.")))
    }
    else { stop(" unknown estimation method"); }
    ptox <- prior^exp(est)
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- prior^exp(ub)
      ptoxU <- prior^exp(lb)
    }
  }
  else if (model=="logistic") {
    dosescaled <- log(prior/(1-prior)) - intcpt
    if (!all(dosescaled<0)) {
      stop( "Intercept parameter in logit model is too small: scaled doses > 0!")
    }

    x1p <- dosescaled[level[include]]

    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrmlgt,c(-10,10),x1p,y1p,w1p,intcpt,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2lgt,-Inf,Inf,x1p,y1p,w1p,500,intcpt,abs.tol=0)[[1]] / integrate(crmhlgt,-Inf,Inf,x1p,y1p,w1p,500,intcpt,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      den <- integrate(crmhlgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]]
      est <- integrate(crmhtlgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den
      if (var.est) { e2 <- integrate(crmht2lgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den; }
    }
    else { stop(" unknown estimation method"); }
    ptox <- (1 + exp(-intcpt-exp(est)*dosescaled))^{-1}
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- (1 + exp(-intcpt-exp(ub)*dosescaled))^{-1}
      ptoxU <- (1 + exp(-intcpt-exp(lb)*dosescaled))^{-1}
    }
  }
  else { stop(" model specified not available."); }

  if (all(ptox<=target)) { rec <- length(prior); }
  else if (all(ptox>=target)) { rec <- 1; }
  else { rec <- order(abs(ptox-target))[1]; }
  if (!var.est) { post.var <- ptoxU <- ptoxL <- NA; }
  foo <- list(prior=prior, target=target, tox=tox, level=level,
              dosename=dosename, subset=pid[include], estimate=est,
              model=model, prior.var=scale^2, post.var=post.var,method=method,
              mtd=rec, include=include, pid=pid, model.detail=model.detail,intcpt=intcpt,
              ptox=ptox, ptoxL=ptoxL, ptoxU=ptoxU, conf.level=conf.level,
              patient.detail=patient.detail,tite=FALSE,dosescaled=dosescaled,var.est=var.est)
  class(foo) <- "mtd"
  foo
}

