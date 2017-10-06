test3 <- function(dim, fun, doe.size, nsteps = 10, criter = 'EI', 
                  xopt = NULL, fopt = NULL,
                  lower, upper, cov.type = "matern3_2", 
                  verbose = FALSE, opts = list()) {
  
  design.fact <- matrix(runif(doe.size * dim), doe.size, dim) * 
    kronecker(matrix(1, doe.size, 1), t(upper - lower)) + 
    kronecker(matrix(1, doe.size, 1), t(lower)) %>%
    as.data.frame %>%
    set_colnames(paste0('x', 1:dim))
  
  response.branin <- apply(design.fact, 1, fun)
  response.branin <- data.frame(response.branin)
  names(response.branin) <- "y"
  
  # model identification
  # fitted.model1 <- km(~1, design = design.fact, response=response.branin,
  #                     nugget = 1e-8,
  #                     covtype = "gauss", control=list(pop.size=50, 
  #                                                     trace=FALSE))
  # # EGO n steps
  # oEGO <- EGO.nsteps(model = fitted.model1, fun=fun, nsteps=nsteps,
  #                    lower = lower, upper=upper, control=list(pop.size=20, 
  #                                                           print.level = 0,
  #                                                           BFGSburnin=2))
  
  model <- km(~1, design = design.fact, response = response.branin, 
                      cov.type, 
              lower = 1e-10 * (upper - lower),
              upper = 1e2 * (upper - lower),
              optim.method = 'BFGS',
              nugget = 1e-8,
              nugget.estim = T,
              multistart = 10,
              control = list(trace = F,
                             maxit = 2e3 * dim))
  
  # run EGO for nsteps
  oEGO <- EGO.nsteps(model, fun, nsteps, lower, upper,
              control = list(trace = FALSE,
                             pop.size = 20, 
                             print.level = 0,
                             max.generations = 200,
                             wait.generations = 10,
                             BFGSburnin = 5))
  
  print(oEGO$par)
  print(oEGO$value)
  
  out <- list()
  # calculate the distance to the optima in search/objective space
  f_dist_best <- rep(0, nsteps)
  x_dist_best <- rep(0, nsteps)
  
  new_y <- oEGO$value
  par <- oEGO$par
  
  for (i in seq(nsteps)) {
    f_dist_best[i] <- min(min(response.branin), new_y[1:i])
  }
  out[['f_dist_best']] <- f_dist_best
  
  if (!is.null(xopt)) {
    if (is.character(xopt)) {
      xopt <- eval(parse(text = xopt))
    }
    
    x_dist <- function(xopt, x) {
      nrow <- nrow(x)
      xopt <- kronecker(matrix(1, nrow, 1), t(xopt))
      rowSums((x - xopt) ^ 2) %>% sqrt
    }
    
    if (is.null(dim(xopt))) {
      dist <- x_dist(par, xopt)
      dist.init <- x_dist(design.fact, xopt) %>% min
    } else {
      dist <- apply(xopt, 1, x_dist, x = par) %>% 
        apply(MARGIN = 1, FUN = min)
      dist.init <- apply(xopt, 1, x_dist, x = design.fact) %>% 
        apply(MARGIN = 1, FUN = min) %>% min
    }
    
    for (i in seq(nsteps)) {
      x_dist_best[i] <- min(c(dist.init, dist[1:i]))
    }
    out[['x_dist_best']] <- x_dist_best
  }
  
  out[['par']] <- par
  out
}

test2 <- function(dim, fun, doe.size, nsteps = 10, criter = 'EI', 
                  xopt = NULL, fopt = NULL,
                  lower, upper, cov.type = "matern3_2", verbose = FALSE, 
                  opts = list()) {
  
  design.fact <- matrix(runif(doe.size * dim), doe.size, dim) * 
    kronecker(matrix(1, doe.size, 1), t(upper - lower)) + 
    kronecker(matrix(1, doe.size, 1), t(lower)) %>%
    as.data.frame %>%
    set_colnames(paste0('x', 1:dim))
  
  response.branin <- apply(design.fact, 1, fun)
  response.branin <- data.frame(response.branin)
  names(response.branin) <- "y"
  
  y.mean <- apply(response.branin, 2, mean)
  y.sd <- apply(response.branin, 2, sd)
  response.branin.norm <- (response.branin - y.mean) / y.sd
  
  # model identification
  # model <- km(~1, design = design.fact, response = response.branin.norm,
  #             covtype = "matern3_2", nugget.estim = T, 
  #             control = list(pop.size = 50, trace = FALSE))
  # # EGO n steps
  # oEGO <- EGO(model = model, fun = fun, nsteps = nsteps,
  #             criter = 'EI',
  #             y.mean = y.mean, y.sd = y.sd,
  #             lower = lower, upper = upper, 
  #             control = list(pop.size = 20, 
  #                            print.level = 0,
  #                            BFGSburnin = 2))
  # 
  model <- km(~1, design = design.fact, response = response.branin.norm, cov.type, 
              lower = 1e-10 * (upper - lower),
              upper = 1e2 * (upper - lower),
              optim.method = 'BFGS',
              nugget = 1e-8,
              nugget.estim = T,
              multistart = 10,
              control = list(trace = F,
                             maxit = 2e3 * dim))
  
  # run EGO for nsteps
  oEGO <- EGO(model, fun, nsteps, lower, upper,
              criter = criter,
              y.mean = y.mean, y.sd = y.sd,
              control = list(trace = FALSE,
                             pop.size = 20, 
                             print.level = 0,
                             max.generations = 200,
                             wait.generations = 10,
                             BFGSburnin = 5),
              verbose = verbose)
  
  print(oEGO$par)
  print(oEGO$value)
  
  out <- list()
  # calculate the distance to the optima in search/objective space
  f_dist_best <- rep(0, nsteps)
  x_dist_best <- rep(0, nsteps)
  
  new_y <- oEGO$value
  par <- oEGO$par
  
  for (i in seq(nsteps)) {
    f_dist_best[i] <- min(min(response.branin), new_y[1:i])
  }
  out[['f_dist_best']] <- f_dist_best
  
  if (!is.null(xopt)) {
    if (is.character(xopt)) {
      xopt <- eval(parse(text = xopt))
    }
    
    x_dist <- function(xopt, x) {
      nrow <- nrow(x)
      xopt <- kronecker(matrix(1, nrow, 1), t(xopt))
      rowSums((x - xopt) ^ 2) %>% sqrt
    }
    
    if (is.null(dim(xopt))) {
      dist <- x_dist(par, xopt)
      dist.init <- x_dist(design.fact, xopt) %>% min
    } else {
      dist <- apply(xopt, 1, x_dist, x = par) %>% 
        apply(MARGIN = 1, FUN = min)
      dist.init <- apply(xopt, 1, x_dist, x = design.fact) %>% 
        apply(MARGIN = 1, FUN = min) %>% min
    }
    
    for (i in seq(nsteps)) {
      x_dist_best[i] <- min(c(dist.init, dist[1:i]))
    }
    out[['x_dist_best']] <- x_dist_best
  }
  
  out[['par']] <- par
  out
}

test <- function(dim, fun, doe.size, nsteps, lower, upper,
                 criter = 'EI', xopt = NULL, fopt = NULL,
                 cov.type = "matern3_2", verbose = FALSE, opts = list()) {
  
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
  
  X <- matrix(runif(doe.size * dim), doe.size, dim) * 
    kronecker(matrix(1, doe.size, 1), t(upper - lower)) + 
    kronecker(matrix(1, doe.size, 1), t(lower)) %>%
    as.data.frame %>%
    set_colnames(paste0('x', 1:dim))
  
  y <- data.frame(apply(X, 1, fun)) %>%
    set_names('y')
  
  # Normalization of the output is needed for MGF
  y.mean <- apply(y, 2, mean)
  y.sd <- apply(y, 2, sd)
  y.norm <- (y - y.mean) / y.sd
  
  # model identification
  model <- km(~1, design = X, response = y.norm, cov.type, 
              lower = 1e-10 * (upper - lower),
              upper = 1e2 * (upper - lower),
              optim.method = 'BFGS',
              nugget.estim = T,
              multistart = 3,
              control = list(trace = F,
                             maxit = 1e2 * dim))
  
  # run EGO for nsteps
  oEGO <- EGO(model, fun, nsteps, lower, upper,
              criter = criter,
              y.mean = y.mean, y.sd = y.sd,
              control = list(trace = FALSE,
                             pop.size = 20, 
                             print.level = 0,
                             max.generations = 200,
                             wait.generations = 3,
                             BFGSburnin = 5),
              verbose = verbose)
  
  print(oEGO$par)
  print(oEGO$value)
  
  out <- list()
  # calculate the distance to the optima in search/objective space
  f_dist_best <- rep(0, nsteps)
  x_dist_best <- rep(0, nsteps)
  
  new_y <- oEGO$value
  par <- oEGO$par
  
  for (i in seq(nsteps)) {
    f_dist_best[i] <- min(min(y), new_y[1:i])
  }
  out[['f_dist_best']] <- f_dist_best
  
  if (!is.null(xopt)) {
    if (is.character(xopt)) {
      xopt <- eval(parse(text = xopt))
    }
    
    x_dist <- function(xopt, x) {
      nrow <- nrow(x)
      xopt <- kronecker(matrix(1, nrow, 1), t(xopt))
      rowSums((x - xopt) ^ 2) %>% sqrt
    }
    
    if (is.null(dim(xopt))) {
      dist <- x_dist(xopt, par)
      dist.init <- x_dist(xopt, X) %>% min
    } else {
      dist <- apply(xopt, 1, x_dist, x = par) %>% 
        apply(MARGIN = 1, FUN = min)
      dist.init <- apply(xopt, 1, x_dist, x = X) %>% 
        apply(MARGIN = 1, FUN = min) %>% min
    }
    
    for (i in seq(nsteps)) {
      x_dist_best[i] <- min(dist.init, dist[1:i])
    }
    out[['x_dist_best']] <- x_dist_best
  }
  out[['par']] <- par
  out
}
