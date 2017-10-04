library(DiceOptim)
library(DiceKriging)
set.seed(123)

d <- 2
n <- 9
fun <- branin

design.fact <- expand.grid(seq(-5,10,length=3), seq(0,15,length=3))
names(design.fact)<-c("x1", "x2")
design.fact <- data.frame(design.fact)
names(design.fact)<-c("x1", "x2")
response.branin <- apply(design.fact, 1, fun)
response.branin <- data.frame(response.branin)
names(response.branin) <- "y"
# model identification
fitted.model1 <- km(~1, design=design.fact, response=response.branin,
                    covtype="gauss", control=list(pop.size=50, 
                                                  trace=FALSE), 
                    parinit=c(0.5, 0.5) * 15)
# EGO n steps
library(rgenoud)
nsteps <- 5 # Was 10, reduced to 5 for speeding up compilation
lower <- c(-5, 0)
upper <- c(10, 15)
oEGO <- EGO.nsteps(model=fitted.model1, fun=fun, nsteps=nsteps,
                   lower=lower, upper=upper, control=list(pop.size=20, 
                                                          print.level = 0,
                                                          BFGSburnin=2))
print(oEGO$par)
print(oEGO$value)

n.grid <- 15 # Was 20, reduced to 15 for speeding up compilation
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, branin)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid, 40)
title("Branin function")
points(design.fact[,1], design.fact[,2], pch=17, col="blue")
points(oEGO$par, pch=19, col="red")
text(oEGO$par[,1], oEGO$par[,2], labels=1:nsteps, pos=3)


