library(DiceKriging)
library(DiceOptim)
library(lhs)
library(ggplot2)
library(grid)
library(gridExtra)
library(magrittr)

# generate the gradient field
# Set the working directory
setwd('~/Dropbox/code_base/acquisition/R')
source('./criteria.R')

set.seed(1234)

dim <- 2
n_sample <- 50

# target function
func <- branin

# Boundary of the function
lower <- c(-5, 0)
upper <- c(10, 15)

# Latin Hyper Cube Sampling
design.lhs <- as.data.frame(optimumLHS(n = n_sample, k = dim) * 
                       kronecker(matrix(1, n_sample, 1), t(upper - lower))
                     + kronecker(matrix(1, n_sample, 1), t(lower))) 

colnames(design.lhs) <- paste0('x', 1:dim)

respones <- apply(design.lhs, 1, func) %>% 
  {(.  - mean(.)) / sqrt(var(.))} %>%
  as.data.frame
names(respones) <- "y"

design.lhs %<>% apply(MARGIN = 2, FUN = . %>% {(.  - mean(.)) / sqrt(var(.))})

# build an Ordinary Kriging model
cov.type <- 'gauss'
gp <- km(~1, design = design.lhs, response = respones, 
         covtype = cov.type, nugget.estim = T, multistart = 10,
         control = list(trace = F,
                        pgtol = 1e-15, 
                        factr = 1e4,
                        maxit = 1e5))

# criterion and its gradient: Moment Generating functions
criter <- . %>% MGF(model = gp, t = 0.5)
criter.dx <- . %>% MGF.dx(model = gp, t = 0.5)
if (1 < 2) {
  criter <- . %>% PI(model = gp)
  criter.dx <- . %>% PI.dx(model = gp)
}

# criter <- . %>% {predict(object = gp, newdata = ., type = "UK", 
#                         checkNames = FALSE, se.compute = TRUE, 
#                         cov.compute = FALSE)$mean}

# criter.dx <- . %>% grad(model = gp) %>% '$'(kriging.sd.grad)

# generate the grid for plotting
step <- 0.1
x <- seq(lower[1], upper[1], by = step)
y <- seq(lower[2], upper[2], by = step)
grid <- expand.grid(x, y)

len.x <- length(x)
len.y <- length(y)

X <- matrix(grid$Var1, len.y, len.x, byrow = F)
Y <- matrix(grid$Var2, len.y, len.x, byrow = F)

# calculate the citerion 
values <- apply(grid, 1, criter) %>%
# values <- predict(object = gp, newdata = grid, type = "UK", 
#                   checkNames = FALSE, se.compute = TRUE, 
#                   cov.compute = FALSE)$sd %>%
  matrix(len.y, len.x, byrow = F)

# p1 <- data.frame(x = grid[, 1], y = grid[, 2], 
#                  criter = matrix(values, len.x * len.y, 1)) %>%
#   ggplot(aes(x, y, z = criter)) + 
#   geom_contour(aes(colour = ..level..), bins = 80) +
#   scale_color_gradient(low = "white", high = "black") +
#   xlim(lower[1], upper[1])
# 
# print(p1)

# and the gradient
step <- 0.5
x <- seq(lower[1], upper[1], by = step)
y <- seq(lower[2], upper[2], by = step)
grid <- expand.grid(x, y)

len.x <- length(x)
len.y <- length(y)

X1 <- matrix(grid$Var1, len.y, len.x, byrow = F)
Y1 <- matrix(grid$Var2, len.y, len.x, byrow = F)

gradient <- apply(grid, 1, criter.dx)

dx <- matrix(gradient[1, ], len.y, len.x, byrow = F)
dy <- matrix(gradient[2, ], len.y, len.x, byrow = F)

save(X, Y, X1, Y1, values, dx, dy, file = './grad.RData')
