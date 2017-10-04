library(DiceOptim)
library(magrittr)

# options(warn = -1) # switch warnings off

setwd('~/Dropbox/code_base/acquisition/R')

wd <- getwd()
source('EGO.R')
source('fitness.R')
source('call_EGO.R')

# ------------------------------------ test settings -----------------------------
dims <- c(2)
ndim <- length(dims)

doe.size <- 21
nsteps <- 10
nrun <- 15

test.mode <- c('OK')
nalgorithm <- length(test.mode)

# fun.list <- c('ackley', 'rastrigin', 'schwefel', 'griewank')
fun.list <- c('branin')
nfun <- length(fun.list)
cov.type <- "matern3_2"

hist_best <- array(NA, dim = c(nrun, nsteps + 1, nalgorithm))

for (j in seq_along(dims)) {
  dim <- dims[j]
  
  for (k in seq_along(fun.list)) {
    
    fun.name <- fun.list[k]
    fun <- get(fun.name)
    
    fun.attr <- attributes(fun)
    fopt <- fun.attr$fopt
    xopt <- fun.attr$xopt
    lower <- eval(parse(text = fun.attr$lower))
    upper <- eval(parse(text = fun.attr$upper))
    
    for (i in seq_along(test.mode)) {
      
      # seeds <- ceiling(1000*abs(rnorm(nrun)))
      
      cat(paste('dim:', dim, 'function:', fun.name, "fopt:", fopt), '\n')
      
      for (k in seq(nrun)) {
        cat('run', k, '\n')
        
        # set.seed(123)
        # res <- test2(dim, fun, doe.size, nsteps, criter = 'EI', xopt = xopt,
        #             fopt = fopt, lower, upper, cov.type, verbose = TRUE)
        # f_dist <- res$f_dist_best
        # x_dist <- res$x_dist_best
        # 
        # cat('f_dist:', f_dist, '\n')
        # cat('x_dist:', x_dist, '\n')
        # 
        # # validation
        # set.seed(123)
        # 
        # res <- test3(dim, fun, doe.size, nsteps, criter = 'EI', xopt = xopt,
        #             fopt = fopt, lower, upper, cov.type, verbose = F)
        # f_dist <- res$f_dist_best
        # x_dist <- res$x_dist_best
        # 
        # cat('f_dist:', f_dist, '\n')
        # cat('x_dist:', x_dist, '\n')
        # browser()
        
        res <- test(dim, fun, doe.size, nsteps, criter = 'MGF', xopt = xopt,
                    fopt = fopt, lower, upper, cov.type, verbose = F)
        f_dist <- res$f_dist_best
        x_dist <- res$x_dist_best
        
        cat('f_dist:', f_dist, '\n')
        cat('x_dist:', x_dist, '\n')
        
        browser()
        hist_best[k, , i] <- as.vector(f_dist)
      }
    }
    
    # data processing and plotting and calculate performance metrics
    SRE <- abs(hist_best - fopt)  # squared relative error
    SRE.mean <- apply(SRE, c(2, 3), mean) %>% t
    SRE.sd <- apply(SRE, c(2, 3), sd) %>% t
    
    cat('MSRE:', SRE.mean, '\n')
    
    #save.time <- format(Sys.time(), '%Y-%m-%d_%H:%M:%S')
    data.path <- file.path(wd, 'data', fun.name)
    print(data.path)
    dir.create(data.path, recursive = TRUE, showWarnings = FALSE)
    
    for (i in seq_along(test.mode)) {
      csv.name <- paste0(dim, "D_", nrun, 'run', '.csv')  
      
      data <- data.frame(steps = 1:nsteps, y = SRE.mean[i, ], se = SRE.sd[i, ],
                         alg.name = rep(test.mode[i], nsteps))
      write.csv(x = data, file = file.path(data.path, csv.name), row.names = FALSE)
    }
  }
}

