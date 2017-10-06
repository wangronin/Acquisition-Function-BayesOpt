library(DiceOptim)
library(magrittr)

# options(warn = -1) # switch warnings off

setwd('~/Dropbox/code_base/acquisition/R')

wd <- getwd()
source('EGO.R')
source('fitness.R')
source('example.R')

# ------------------------------------ test settings -----------------------------
dims <- c(2)
ndim <- length(dims)

set.seed(123)

doe.size <- 50
nsteps <- 30
nrun <- 30

test.mode <- c('PI')
nalgorithm <- length(test.mode)

# fun.list <- c('ackley', 'rastrigin', 'schwefel', 'griewank')
fun.list <- c('rastrigin')
nfun <- length(fun.list)
cov.type <- "matern3_2"

hist_best_f <- array(NA, dim = c(nrun, nsteps, nalgorithm))

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
      
      for (n in seq(nrun)) {
        cat('run', n, '\n')
        
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
        
        res <- test(dim, fun, doe.size, nsteps, criter = test.mode[i],
                    xopt = xopt, fopt = fopt, lower, upper, cov.type, verbose = F)
        f_dist <- res$f_dist_best
        x_dist <- res$x_dist_best
        
        # cat('f_dist:', f_dist, '\n')
        # cat('x_dist:', x_dist, '\n')
        hist_best_f[n, , i] <- f_dist
      }
    }
    
    # data processing and plotting and calculate performance metrics
    SRE <- abs(hist_best_f - fopt)  # squared relative error
    SRE.mean <- apply(SRE, 2, mean) %>% t
    SRE.sd <- apply(SRE, 2, sd) %>% t
    
    data.path <- file.path(wd, 'data', fun.name)
    dir.create(data.path, recursive = TRUE, showWarnings = FALSE)
    
    for (i in seq_along(test.mode)) {
      criter <- test.mode[i]
      cat(criter, '\n')
      cat('MSRE:', SRE.mean[i, ], '\n')
      
      data <- data.frame(t(hist_best_f[, , i])) %>% 
        `colnames<-`(paste0('run', seq(nrun)))
      
      csv.name <- paste0(criter, '_', dim, "D", '.csv')
      write.csv(data, file = file.path(data.path, csv.name), row.names = FALSE)
    }
  }
}

