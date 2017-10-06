eps <- 0

# TODO: clean the code and 
rastrigin <- function(xx)
{
  d <- length(xx)
  
  sum <- sum(xx^2 - 10*cos(2*pi*xx))
  
  y <- 10*d + sum + eps
  return(y)
}

attributes(rastrigin) <- list(fopt = eps, 
                              xopt = 'rep(0, dim)',
                              lower = 'rep(-5.12, dim)', 
                              upper = 'rep(5.12, dim)')

schwefel <- function(xx)
{
  d <- length(xx)
  sum <- sum(xx*sin(sqrt(abs(xx))))
  y <- 418.9828872724339 * d - sum + eps
  
  return(y)
}

attributes(schwefel) <- list(fopt = eps, 
                             xopt = 'rep(420.96874636, dim)',
                             lower = 'rep(-500, dim)', 
                             upper = 'rep(500, dim)')

griewank <- function(xx)
{
  ii <- c(1:length(xx))
  sum <- sum(xx^2/4000)
  prod <- prod(cos(xx/sqrt(ii)))
  
  y <- sum - prod + 1 + eps
  return(y)
  
}

attributes(griewank) <- list(fopt = eps,
                             xopt = 'rep(0, dim)',
                             lower = 'rep(-600, dim)', 
                             upper = 'rep(600, dim)')

ackley <- function(xx, a = 20, b = 0.2, c = 2*pi){
  
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  y <- term1 + term2 + a + exp(1) + eps
  return(y)
}

attributes(ackley) <- list(xopt = 'rep(0, dim)', 
                           fopt = eps, 
                           lower = 'rep(-15, dim)', 
                           upper = 'rep(30, dim)')

himmelblau <- function(xx){
	x1 <- xx[1]
	x2 <- xx[2]
	(x1 ^ 2 + x2 - 11) ^ 2 +  (x1 + x2 ^ 2 - 7) ^ 2 + eps
} 

attributes(himmelblau) <- list(xopt = rbind(c(3., 2.), 
                                            c(-2.805118, 3.131312),
                                            c(-3.779310, -3.283186), 
                                            c(3.584428, -1.848126)), 
                               fopt = eps, 
                               lower = c(-5, -5), 
                               upper = c(5, 5))

rosenbrock <- function(xx)
{
  
  d <- length(xx)
  xi <- xx[1:(d-1)]
  xnext <- xx[2:d]
  
  sum <- sum(100 * (xnext - xi ^ 2) ^ 2 + (xi - 1) ^ 2)
  y <- sum
  
  return(y)
}

attributes(rosenbrock) <- list(xopt = 'rep(1, dim)', 
                               fopt = eps, 
                               lower = 'rep(-5, dim)', 
                               upper = 'rep(5, dim)')

# TODO: check the xopt, fopt for those functions
branin <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 -
    1/(8 * pi)) * cos(x1) + 10
}

attributes(branin) <- list(fopt = 0.397887,
                           xopt = rbind(c(-pi, 12.275), c(pi, 2.275),
                                        c(9.42478, 2.475)),
                           lower = 'c(-5, 0)',
                           upper = 'c(10, 15)')

# branin function rescaled to the unit box
branin2 <- function(x) {
  x1 <- x[1] * 15 - 5
  x2 <- x[2] * 15
  (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 
     1/(8 * pi)) * cos(x1) + 10
  
}

attributes(branin2) <- list(fopt = 0.397887,
                           xopt = ((rbind(c(-pi, 12.275), c(pi, 2.275),
                                        c(9.42478, 2.475)) - c(-5, 0))) / 15,
                           lower = 'c(0, 1)', 
                           upper = 'c(0, 1)')

M <- function(xx, alpha = 6){
  n <- length(xx)
  - sum(sin(5*pi*xx)^alpha) / n + 11
}


diffpow <- function(x){
  d <- length(x)
  x <- abs(x)
  p <- (10* seq(0, d-1) / (d-1) + 2)
  f <- sum(x ^ p)
  return(f)
}
