# implement those as classes
# compute the auxiliary variables 
aux <- function(x, model, plugin, type, minimization) {
  if (is.null(plugin)) {
    if (minimization) {
      plugin <- min(model@y)
    }
    else {
      plugin <- -max(model@y)
    }
  }
  
  ymin <- plugin
  d <- length(x)
  
  if (d != model@d) stop("x does not have the right size")
  
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  predx <- predict(object = model, newdata = newdata, type = type, 
                   checkNames = FALSE)
  
  return(list(predx$mean, predx$sd, ymin))
}

# compute the gradient of Kriging mean and std. w.r.t the input vector
km.dx <- function(x, model, type = "UK", envir = NULL) {
  
  d <- length(x)
  if (d != model@d)
    stop("x does not have the right size")
  
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  
  T <- model@T
  X <- model@X
  z <- model@z
  u <- model@M
  covStruct <- model@covariance
  
  # calculate the required variables
  if (is.null(envir)) {
    predx <- predict(object = model, newdata = newdata, type = type, 
                     checkNames = FALSE, se.compute = TRUE, cov.compute = FALSE)
    kriging.sd <- predx$sd
    v <- predx$Tinv.c
    c <- predx$c
  } else { # no prediction needed if the enviroment is provided
    toget <- matrix(c("kriging.sd", "c", "Tinv.c"), 1, 3)
    apply(toget, 2, get, envir = envir)
    kriging.sd <- envir$kriging.sd
    c <- envir$c
    v <- envir$Tinv.c
  }
  
  F.newdata <- model.matrix(model@trend.formula, data = newdata)
  dc <- covVector.dx(x = newdata.num, X = X, object = covStruct, 
                     c = c)
  
  # test for covVector gradient function
  # dc.test <- covVector.dx.test(x = newdata.num, X = X, object = covStruct, 
  # c = c)
  f.deltax <- trend.deltax(x = newdata.num, model = model)
  W <- backsolve(t(T), dc, upper.tri = FALSE)
  kriging.mean.grad <- t(W) %*% z + t(model@trend.coef %*% f.deltax)
  
  if (type == "UK") {
    tuuinv <- solve(t(u) %*% u)
    kriging.sd2.grad <- t(-2 * t(v) %*% W +
                            2 * (F.newdata - t(v) %*% u) %*% tuuinv %*%
                            (f.deltax - t(t(W) %*% u)))
    
  } else {
    kriging.sd2.grad <- t(-2 * t(v) %*% W)
  }
  
  kriging.sd.grad <- kriging.sd2.grad / (2 * kriging.sd)
  
  return(list(kriging.mean.grad = kriging.mean.grad, 
              kriging.sd.grad = kriging.sd.grad))
}

# The moment generating function of the improvement
# Note that this function can generate infinite values
MGF <- function(x, model, plugin = NULL, t = 1, type = "UK",
                minimization = TRUE, envir = NULL) {
  
  # tmp <- aux(x, model, plugin, type, minimization)
  if (is.null(plugin)) {
    if (minimization) {
      plugin <- min(model@y)
    }
    else {
      plugin <- -max(model@y)
    }
  }
  
  d <- length(x)
  
  if (d != model@d) {
    stop("x does not have the right size")
  }
  
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  predx <- predict(object = model, newdata = newdata, type = type, 
                   checkNames = FALSE)
  
  m <- kriging.mean <- predx$mean
  s <- kriging.sd <- predx$sd
  ymin <- plugin
  
  if (s/sqrt(model@covariance@sd2) < 1e-06) 
    return(0)
  
  m.prime <- m - t * s ^ 2
  
  beta <- (ymin - m) / s
  beta.prime <- (ymin - m.prime) / s
  
  if (!minimization) {
    m <- -m
  }
  
  # Note that although the MGF-based acquisition function always exists
  # for finite t theoretically, the exponential term below could explode
  # in practice. Therefore, the normlization of y is required
  term <- t * (ymin - m - 1)
  res <- pnorm(beta.prime) * exp(term + t ^ 2 * s ^ 2 / 2)
  
  # for debug purposes: numerical explosion
  if (is.infinite(res)) 
    browser()
  
  if (!is.null(envir)) {
    assign('kriging.mean', kriging.mean, envir = envir)
    assign("kriging.sd", kriging.sd, envir = envir)
    assign("c", predx$c, envir = envir)
    assign("Tinv.c", predx$Tinv.c, envir = envir)
  }
  res
}

# Compute the gradient of the moment generating function based acquisition function
MGF.dx <- function(x, model, plugin = NULL, t = 1, type = "UK", 
                   minimization = TRUE, envir = NULL) {
  
  if (is.null(plugin)) {
    if (minimization) {
      plugin <- min(model@y)
    } 
    else {
      plugin <- -max(model@y)
    }
  }
  
  m <- plugin
  d <- length(x)
  if (d != model@d)
    stop("x does not have the right size")
  
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) <- colnames(model@X)
  ymin <- plugin
  
  if (is.null(envir)) {
    predx <- predict(object = model, newdata = newdata, type = type, 
                     checkNames = FALSE, se.compute = TRUE, cov.compute = FALSE)
    
    kriging.mean <- predx$mean
    if (!minimization) 
      kriging.mean <- -kriging.mean
    
    kriging.sd <- predx$sd
    Tinv.c <- predx$Tinv.c
    c <- predx$c
    
  } else {
    toget <- matrix(c("kriging.sd", "kriging.mean", "c", "Tinv.c"), 1, 4)
    apply(toget, 2, get, envir = envir)
    kriging.mean <- envir$kriging.mean
    kriging.sd <- envir$kriging.sd
    c <- envir$c
    Tinv.c <- envir$Tinv.c
  }
  
  m <- kriging.mean
  s <- kriging.sd
  m.prime <- m - t * s ^ 2
  
  beta <- (ymin - m) / s
  beta.prime <- (ymin - m.prime) / s
  
  if (s/sqrt(model@covariance@sd2) < 1e-06) {
    dx <- rep(0, d)
  } else {
    # Compute the y^hat and sd gradient and attach to the running frame
    list2env(km.dx(x, model, type = type, envir = environment()), 
             envir = environment())
    
    if (!minimization) 
      kriging.mean.grad <- -kriging.mean.grad
    
    m.dx <- kriging.mean.grad
    s.dx <- kriging.sd.grad
    
    term1 <- exp(ymin * t + t ^ 2 * s ^ 2 / 2 - m * t - t)
    m.prime.dx <- m.dx - 2 * t * s * s.dx
    beta.prime.dx <- -(m.prime.dx + beta.prime * s.dx) / s
    
    dx <- dnorm(beta.prime) * term1 * beta.prime.dx + 
      pnorm(beta.prime) * term1 * ((t ^ 2) * s * s.dx - t * m.dx)
    
  }
  dx
}

# The Generalized Expected Improvement
# The n-th order moment of the improvement, about the origin
GEI <- function(x, model, plugin = NULL, g = 2, type = "UK", minimization = TRUE, 
                envir = NUL) {
  
  tmp <- aux(x, model, plugin, type, minimization)
  
  y.hat <- tmp[[1]]
  s <- tmp[[2]]
  ymin <- tmp[[3]]
  
  u <- (ymin - y.hat) / s
  t0 <-  pnorm(u)
  t1 <- -dnorm(u)
  
  if (s/sqrt(model@covariance@sd2) < 1e-06) return(0)
  
  if (g == 0) return(t0)
  else if (g == 1) return(s*(u * t0 - t1))
  else {
    t <- c(t0, t1, rep(0, g - 1))
    for (k in seq(3, g + 1)) {
      t[k] <- u ^ (k - 2) * t1 + (k - 2) * t[k - 2]
    }
    
    k <- seq(0, g)
    res <- (s ^ g) * sum(((-1) ^ k) * choose(g, k) * (u ^ (g - k)) * t)
  }
  
  if (res < 0) res <- 0 # correct slightly negative values
  
  return(res)
}

PI <- function(x, model, plugin = NULL, type = "UK", minimization = TRUE, 
               envir = NULL) {
  
  if (is.null(plugin)) {
    if (minimization) {
      plugin <- min(model@y)
    }
    else {
      plugin <- -max(model@y)
    }
  }
  
  m <- plugin
  d <- length(x)
  
  if (d != model@d) {
    stop("x does not have the right size")
  }
  
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  predx <- predict(object = model, newdata = newdata, type = type, 
                   checkNames = FALSE)
  
  kriging.mean <- predx$mean
  
  if (!minimization) {
    kriging.mean <- -kriging.mean
  }
  
  kriging.sd <- predx$sd
  xcr <- (m - kriging.mean)/kriging.sd
  
  if (kriging.sd/sqrt(model@covariance@sd2) < 1e-06) {
    res <- 0
    xcr <- xcr.prob <- NULL
  }
  else {
    res <- xcr.prob <- pnorm(xcr)
  }
  
  if (!is.null(envir)) {
    assign("xcr", xcr, envir = envir)
    assign("xcr.prob", xcr.prob, envir = envir)
    assign('kriging.mean', kriging.mean, envir = envir)
    assign("kriging.sd", kriging.sd, envir = envir)
    assign("c", predx$c, envir = envir)
    assign("Tinv.c", predx$Tinv.c, envir = envir)
  }
  return(res)
}


PI.dx <- function(x, model, plugin = NULL, type = "UK", minimization = TRUE, 
                  envir = NULL) {
  
  if (is.null(plugin)) {
    if (minimization) {
      plugin <- min(model@y)
    } 
    else {
      plugin <- -max(model@y)
    }
  }
  
  m <- plugin
  d <- length(x)
  if (d != model@d)
    stop("x does not have the right size")
  
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) <- colnames(model@X)
  
  if (is.null(envir)) {
    predx <- predict(object = model, newdata = newdata, type = type, 
                     checkNames = FALSE, se.compute = TRUE, cov.compute = FALSE)
    kriging.mean <- predx$mean
    if (!minimization) 
      kriging.mean <- -kriging.mean
    
    kriging.sd <- predx$sd
    Tinv.c <- predx$Tinv.c
    c <- predx$c
    xcr <- (m - kriging.mean)/kriging.sd
    
  } else {
    toget <- matrix(c("xcr", "kriging.sd", "kriging.mean", "c", "Tinv.c"), 1, 7)
    apply(toget, 2, get, envir = envir)
    xcr <- envir$xcr
    kriging.mean <- envir$kriging.mean
    kriging.sd <- envir$kriging.sd
    c <- envir$c
    Tinv.c <- envir$Tinv.c
  }
  
  xcr.dens <- if (is.null(xcr)) NULL else dnorm(xcr)
  F.newdata <- model.matrix(model@trend.formula, data = newdata)
  
  if (kriging.sd/sqrt(model@covariance@sd2) < 1e-06) {
    pi.grad <- rep(0, d)
  } else {
    
    # Compute the y^hat and sd gradient and attach to the running frame
    list2env(grad(x, model, type = type, envir = environment()), 
             envir = environment())
    
    if (!minimization) 
      kriging.mean.grad <- -kriging.mean.grad
    
    pi.grad <- -(xcr.dens / kriging.sd) * (kriging.mean.grad 
                                           + xcr * kriging.sd.grad)
  }
  return(pi.grad)
}
