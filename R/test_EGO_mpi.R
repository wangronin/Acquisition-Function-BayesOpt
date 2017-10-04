# EGO testing script using pdbMPI, better than Rmpi: why? I forgot...
suppressMessages(library(pbdMPI))

# initialize MPI 
init()

suppressMessages(library(DiceOptim))
suppressMessages(library(magrittr))
options(warn = -1) # switch warnings off

# get info from the MPI parallel processing
n.run <- comm.size()
rank <- comm.rank()

setwd('~/criteria/')
wd <- getwd()

source('./R/EGO.R')
source('./R/fitness.R')

# ------------------------------------ test settings -----------------------------
dims <- c(2)
n.dim <- length(dims)

doe.size <- 21
n.step <- 100
cov.type <- "matern3_2"

# algorithm variants to be compared
algorithms <- c('EI', 'PI', 'MGF')
n.alg <- length(algorithms)

# test functions
fun.list <- c('branin', 'ackley', 'rastrigin', 'schwefel', 'griewank')
n.fun <- length(fun.list)

verbose <- if (rank == 0) T else F
hist_best <- array(NA, dim = c(n.run, n.step + 1, n.alg))

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
    
    if (rank == 0) {
      cat(paste('test on', dim, 'D', fun.name, 'function, fopt:', fopt, '\n'))
    }
    
    for (i in seq_along(algorithms)) {
      
      criter <- algorithms[i] 
      seeds <- ceiling(1000*abs(rnorm(n.run)))
      
      if (rank == 0) cat(criter, ' ')
      
      # execute the test 
      res <- test(dim, fun, doe.size, n.step, criter = criter, xopt, 
                  lower, upper, cov.type, verbose = verbose)$y_best_hist
      
      # Synchronization
      barrier()
      
      # gather() would collect data from all the nodes in a list
      res.best <- gather(res, rank.dest = 0)
      
      if (rank == 0) {
        
        res.best <- matrix(unlist(res.best), n.run, n.step + 1, byrow = T)
        
        if (nrow(res.best) != n.run) {
          res.best <- t(res.best)
        }
        hist_best[, , i] <- res.best
      }
    }
    
    # On rank 0 node: data processing and plotting
    if (rank == 0) {
      
      # calculate performance metrics
      SRE <- ((hist_best - fopt) / fopt) ^ 2  # squared relative error
      SRE.mean <- apply(SRE, c(2, 3), mean) %>% t
      SRE.sd <- apply(SRE, c(2, 3), sd) %>% t
      
      #save.time <- format(Sys.time(), '%Y-%m-%d_%H:%M:%S')
      data.path <- file.path(wd, 'data', fun.name)
      dir.create(data.path, recursive = TRUE, showWarnings = FALSE)
      
      for (i in seq_along(algorithms)) {
        criter <- algorithms[i]
        cat(criter, '\n')
        cat('MSRE:', SRE.mean[i, ], '\n')
        
        data <- data.frame(t(hist_best[, , i])) %>% 
          `colnames<-`(paste0('run', seq(n.run)))
        
        csv.name <- paste0(criter, '_', dim, "D", '.csv')
        write.csv(data, file = file.path(data.path, csv.name), row.names = FALSE)
      }
      cat('\n\n')
    }
  }
}

finalize()

