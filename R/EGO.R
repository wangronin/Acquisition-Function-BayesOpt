suppressMessages(library(DiceOptim))
suppressMessages(library(rgenoud))
suppressMessages(library(magrittr))
source('criteria.R')

`EGO` <- function(model, fun, nsteps, lower, upper, criterion = 'EI', 
                   parinit = NULL, control = NULL, 
                   kmcontrol = NULL, X.mean = NULL, X.sd = NULL, 
                   y.mean = NULL, y.sd = NULL, verbose = FALSE) {
  
  n <- nrow(model@X)
  if (is.null(kmcontrol$penalty)) 
    kmcontrol$penalty <- model@penalty
  if (length(model@penalty == 0)) 
    kmcontrol$penalty <- NULL
  if (is.null(kmcontrol$optim.method)) 
    kmcontrol$optim.method <- model@optim.method
  if (is.null(kmcontrol$parinit)) 
    kmcontrol$parinit <- model@parinit
  if (is.null(kmcontrol$control)) 
    kmcontrol$control <- model@control
  
  # set initial tempurature...
  if (criterion == 'MGF') 
    t <- 6    
  
  for (i in 1:nsteps) {
    if (verbose) 
      cat(paste0(i, '...'))
    
    criter <- switch(criterion, EI = EI, PI = PI, MGF = MGF)
    criter.dx <- switch(criterion, EI = EI.grad, PI = PI.dx,
                        MGF = MGF.dx)
    
    # maximize the acquisition function
    # sink('/dev/null')
    oEGO <- max_criter(model, criter, criter.dx, lower = lower, 
                       upper = upper, parinit = parinit, control = control)
    # sink(NULL)
    
    # de-normalization and normalization again...
    par <- oEGO$par
    y <- model@y * y.sd + y.mean
    
    # append the new design site to the data set
    model@X <- rbind(model@X, par)
    y <- rbind(y, fun(t(par)))
    y.mean <- apply(y, 2, mean)
    y.sd <- apply(y, 2, sd)
    model@y <- (y - mean(y)) / y.sd
    
    # model hyper-parameters re-estimation
    kmcontrol$parinit <- covparam2vect(model@covariance)
    kmcontrol$control$trace = FALSE
    
    if (model@param.estim) {
      model <- km(formula = model@trend.formula, design = model@X, 
                  response = model@y, covtype = model@covariance@name, 
                  lower = model@lower, upper = model@upper, nugget = NULL, 
                  penalty = kmcontrol$penalty, optim.method = kmcontrol$optim.method, 
                  parinit = kmcontrol$parinit, control = kmcontrol$control, 
                  gr = model@gr, iso = is(model@covariance, "covIso"))
    } else {
      coef.cov <- covparam2vect(model@covariance)
      model <- km(formula = model@trend.formula, design = model@X, 
                  response = model@y, covtype = model@covariance@name, 
                  coef.trend = model@trend.coef, coef.cov = coef.cov, 
                  coef.var = model@covariance@sd2, nugget = NULL, 
                  iso = is(model@covariance, "covIso"))
    }
    
    # decrease the tempurature for MGF and GEI
    if (criterion == 'MGF') {
      t <- t * 0.85   # exponetial decay of the temperature
    }
    
  }
  if (verbose) cat('done\n')
  
  y <- model@y * y.sd + y.mean
  return(list(par = model@X[(n + 1):(n + nsteps), , drop = FALSE], 
              value = y[(n + 1):(n + nsteps), , drop = FALSE], 
              npoints = 1, nsteps = nsteps, lastmodel = model))
}

max_criter <- function(model, criter, criter.dx = NULL, 
                       plugin = NULL, type = "UK", 
                       lower, upper, parinit = NULL, 
                       minimization = TRUE, control = NULL) {
  if (is.null(plugin)) {
        plugin <- min(model@y)
  }
  
  criter.envir <- new.env()
  environment(criter) <- environment(criter.dx) <- criter.envir 
  gr = criter.dx
  d <- ncol(model@X)

  if (is.null(control$print.level)) 
      control$print.level <- 1
  if (d <= 6) 
      N <- 3 * 2^d
  else N <- 32 * d
  if (is.null(control$BFGSmaxit)) 
      control$BFGSmaxit <- N
  if (is.null(control$pop.size)) 
      control$pop.size <- N
  if (is.null(control$solution.tolerance)) 
      control$solution.tolerance <- 1e-21
  if (is.null(control$max.generations)) 
      control$max.generations <- 12
  if (is.null(control$wait.generations)) 
      control$wait.generations <- 2
  if (is.null(control$BFGSburnin)) 
      control$BFGSburnin <- 2
  if (is.null(parinit)) 
      parinit <- lower + runif(d) * (upper - lower)
  domaine <- cbind(lower, upper)

  o <- genoud(criter, nvars = d, max = TRUE, pop.size = control$pop.size, 
              max.generations = control$max.generations, 
              wait.generations = control$wait.generations, 
              hard.generation.limit = TRUE, starting.values = parinit, 
              MemoryMatrix = TRUE, Domains = domaine, default.domains = 10, 
              solution.tolerance = control$solution.tolerance, gr = criter.dx, 
              boundary.enforcement = 2, lexical = FALSE, gradient.check = FALSE, 
              BFGS = TRUE, data.type.int = FALSE, hessian = FALSE, 
              unif.seed = floor(runif(1, max = 10000)), 
              int.seed = floor(runif(1, max = 10000)), 
              print.level = control$print.level, 
              share.type = 0, instance.number = 0, output.path = "stdout", 
              output.append = FALSE, project.path = NULL, P1 = 50, 
              P2 = 50, P3 = 50, P4 = 50, P5 = 50, P6 = 50, P7 = 50, 
              P8 = 50, P9 = 0, P9mix = NULL, BFGSburnin = control$BFGSburnin, 
              BFGSfn = NULL, BFGShelp = NULL, 
              control = list(maxit = control$BFGSmaxit), 
              cluster = FALSE, balance = FALSE, debug = FALSE, model = model, 
              plugin = plugin, type = type, minimization = minimization, 
              envir = criter.envir)
  
  o$par <- t(as.matrix(o$par))
  colnames(o$par) <- colnames(model@X)
  o$value <- as.matrix(o$value)
  colnames(o$value) <- "criterion"
  return(list(par = o$par, value = o$value))
}
