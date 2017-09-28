eps <- 0

M <- function(xx, alpha=6){
	n <- length(xx)
	- sum(sin(5*pi*xx)^alpha) / n + 11
}

rastrigin <- function(xx)
{
  d <- length(xx)
  
  sum <- sum(xx^2 - 10*cos(2*pi*xx))
  
  y <- 10*d + sum + eps
  return(y)
}

attributes(rastrigin) <- list(fopt = eps, 
                              lower = 'rep(-5.12, dim)', 
                              upper = 'rep(5.12, dim)')

schwefel <- function(xx)
{
  d <- length(xx)
  
  sum <- sum(xx*sin(sqrt(abs(xx))))
  
  y <- 418.9829*d - sum + eps
  return(y)
}

attributes(schwefel) <- list(fopt = eps, lower = 'rep(-500, dim)', 
                             upper = 'rep(500, dim)')

griewank <- function(xx)
{
  ii <- c(1:length(xx))
  sum <- sum(xx^2/4000)
  prod <- prod(cos(xx/sqrt(ii)))
  
  y <- sum - prod + 1 + eps
  return(y)
  
}

attributes(griewank) <- list(fopt = eps, lower = 'rep(-600, dim)', 
                             upper = 'rep(600, dim)')

branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi))
{
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1 - t)*cos(x1)
  
  y <- term1 + term2 + s 
  return(y)
}

attributes(branin) <- list(fopt = 0.397887, lower = 'c(-5, 0)', 
                           upper = 'c(10, 15)')

ackley <- function(xx, a = 20, b = 0.2, c = 2*pi){
  
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  y <- term1 + term2 + a + exp(1) + eps
  return(y)
}

attributes(ackley) <- list(xopt = 'rep(0, dim)', fopt = eps, 
                           lower = 'rep(-15, dim)', 
                           upper = 'rep(30, dim)')

himmelblau <- function(xx){
	x1 <- xx[1]
	x2 <- xx[2]
	
	(x1 ^ 2 + x2 - 11) ^ 2 +  (x1 + x2 ^ 2 - 7) ^ 2 + 10
} 

attributes(himmelblau) <- list(xopt='c(3, 2)', fopt=10, lower='c(-5, -5)', upper='c(5, 5)')


rosenbrock <- function(xx)
{
  
  d <- length(xx)
  xi <- xx[1:(d-1)]
  xnext <- xx[2:d]
  
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
	
  y <- sum
  return(y)
}

diffpow <- function(x){
  d <- length(x)
  x <- abs(x)
  p <- (10* seq(0, d-1) / (d-1) + 2)
  f <- sum(x ^ p)
  return(f)
}
