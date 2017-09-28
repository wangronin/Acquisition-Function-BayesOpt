suppressMessages(library(DiceOptim))
suppressMessages(library(rgenoud))
suppressMessages(library(magrittr))

source('./R/criteria.R')

test <- function(dim, fun, doe.size, nsteps = 10, criter = 'EI', lower, upper, 
                 cov.type = "matern5_2", verbose = FALSE, opts = list()) {
  
  # divert the output
  # sink(paste0(Sys.info()['nodename'], '.', Sys.getpid(), '.log'), append = TRUE)
  
  # set the random seed first
  # if ('seed' %in% names(opts)){
  #   set.seed(opts$seed)
  #   set.seed(seed)
  #
  # } else {
  #   #set.seed(as.integer(1000*rank + Sys.getpid() + as.numeric(Sys.time())) %% 1000)
  # }
  
  hist_best <- rep(0, nsteps)
  X <- matrix(runif(doe.size*dim), doe.size, dim) * 
    kronecker(matrix(1, doe.size, 1), t(upper - lower)) + 
    kronecker(matrix(1, doe.size, 1), t(lower)) %>%
    as.data.frame %>%
    set_colnames(paste0('x', 1:dim))
  
  y <- data.frame(apply(X, 1, fun)) %>%
    set_names('y')
  
  # Normalization is needed here for the MGF
  X.mean <- apply(X, 2, mean)
  X.sd <- apply(X, 2, sd)
  y.mean <- apply(y, 2, mean)
  y.sd <- apply(y, 2, sd)
  
  # normalized input and output
  X.norm <- (X - X.mean) / X.sd
  y.norm <- (y - y.mean) / y.sd
  
  lower.norm <- (lower - X.mean) / X.sd
  upper.norm <- (upper - X.mean) / X.sd
  range <- upper.norm - lower.norm
  
  # model identification
  model <- km(~1, X.norm, y.norm, cov.type, 
              lower = 1e-10 * range,
              upper = 1e2 * range,
              optim.method = 'BFGS',
              nugget.estim = TRUE,
              multistart = 20,
              control = list(trace = F,
                             pgtol = 1e-8, 
                             factr = 1e7,
                             maxit = 2e3 * dim))
  
  # run EGO for nsteps
  oEGO <- EGO(model, fun, nsteps, lower, upper,
              criter = criter,
              X.mean = X.mean, X.sd = X.sd,
              y.mean = y.mean, y.sd = y.sd,
              control = list(trace = FALSE,
                             pop.size = 50, 
                             max.generations = 200,
                             wait.generations = 10,
                             BFGSburnin = 5),
              verbose = verbose)
  
  
  new_y <- oEGO$value
  for (i in seq(nsteps)) {
    hist_best[i] <- min(c(t(y), new_y[1:i]))
  }
  return(list(y_best_hist = c(min(y), hist_best), new_X = oEGO$par))
}

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
  
  if (criterion == 'MGF') t <- 7    # tempurature...
  if (verbose) cat('step:')
  
  for (i in 1:nsteps) {
    if (verbose) cat(paste0(i, '. '))
    
    criter <- switch(criterion, EI = EI, PI = PI, MGF = MGF)
    criter.dx <- switch(criterion, EI = EI.grad, PI = PI.dx, MGF = MGF.dx)
    
    # transform the bounds first 
    lower.norm <- (lower - X.mean) / X.sd
    upper.norm <- (upper - X.mean) / X.sd
    
    # maximize the acquisition function
    sink('/dev/null')
    oEGO <- max_criter(model, criter, criter.dx, lower = lower.norm, 
                       upper = upper.norm, parinit = parinit, control = control)
    sink(NULL)
    
    # de-normalization and normalization again...
    par <- oEGO$par * X.sd + X.mean
    X <- model@X * X.sd + X.mean
    y <- model@y * y.sd + y.mean
    
    X <- rbind(X, par)
    y <- rbind(y, fun(t(par)))
    
    X.mean <- apply(X, 2, mean)
    X.sd <- apply(X, 2, sd)
    y.mean <- apply(y, 2, mean)
    y.sd <- apply(y, 2, sd)
    
    # append the new design site to the data set
    model@X <- (X - X.mean) / X.sd
    model@y <- (y - y.mean) / y.sd
    
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
    
    # cooling down for MGF and GEI
    if (criterion == 'MGF') {
      t <- t * 0.85   # exponetial decay of the temperature
    }
  }
  if (verbose) cat('done\n')
  
  return(list(par = X[(n + 1):(n + nsteps), , drop = FALSE], 
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
  f <- criter
  gr <- criter.dx
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
  o <- genoud(f, nvars = d, max = TRUE, pop.size = control$pop.size, 
              max.generations = control$max.generations, 
              wait.generations = control$wait.generations, 
              hard.generation.limit = TRUE, starting.values = parinit, 
              MemoryMatrix = TRUE, Domains = domaine, default.domains = 10, 
              solution.tolerance = control$solution.tolerance, gr = gr, 
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
