### Tutorial for Gaussian Processes, based on
### Roberts et al. (2012). "GPs for time-series modelling"
### This version by Alan Aw (alanaw1), dated: 3/11/2019

## You need these packages 
library(MASS)
library(dplyr)
library(ggplot2)
library(reshape2)

## Build covariance kernel
## Inputs: X1, X2 data,
##        and specified kernel
## Output: matrix of dimension dim(X1) x dim(X2) 
calcK <- function(X1, X2, kernel = "SE", 
                  sigma, h, lambda, alpha, 
                  nu) {
  K <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  # if kernel is WN, specify sigma
  if (kernel == "WN") {
    for (i in 1:nrow(K)) {
      for (j in 1:ncol(K)) {
        K[i,j] <- sigma^2 * as.numeric((i == j))
      }
    }  
  }
  # if kernel is BM (Brownian Motion)
  if (kernel == "BM") {
    for (i in 1:nrow(K)) {
      for (j in 1:ncol(K)) {
        K[i,j] <- min(X1[i],X2[j])
      }
    } 
  }
  # if kernel is BB (Brownian Bridge)
  if (kernel == "BB") {
    for (i in 1:nrow(K)) {
      for (j in 1:ncol(K)) {
        K[i,j] <- min(X1[i],X2[j]) - X1[i]*X2[j]
      }
    } 
  }
  # if kernel is SE, specify h and lambda
  if (kernel == "SE") {
    for (i in 1:nrow(K)) {
      for (j in 1:ncol(K)) {
        K[i,j] <- h^2 * exp(-(X1[i] - X2[j])^2 / lambda^2)
      }
    }
  }
  # if kernel is RQ, specify h, lambda and alpha
  if (kernel == "RQ") {
    for (i in 1:nrow(K)) {
      for (j in 1:ncol(K)) {
        K[i,j] <- h^2 * (1 + (X1[i] - X2[j])^2 / (alpha * lambda^2))^(-alpha)
      }
    }    
  }
  # if kernel is Matern, specify h, lambda and nu
  if (kernel == "Matern") {
    for (i in 1:nrow(K)) {
      for (j in 1:ncol(K)) {
        K[i,j] <- h^2 * 1 / (gamma(nu)*2^(nu - 1)) * 2 * sqrt(nu) * 
          abs(X1[i] - X2[j]) / lambda * Bessel(2 * sqrt(nu) * abs(X1[i] - X2[j]) / lambda, nu)
      }
    }    
  }
  # return K with filled entries
  return(K)
}

## Specify mean function
## Inputs: X is a vector of independent variables
## Output: vector of length dim(X) with means computed 
setMean <- function(X, func = "zero", 
                    lambda, h, 
                    a, b, c, ...) {
  # if func is zero
  if (func == "zero") {
    return(rep(0,length(X)))    
  }
  # if func is exponential decay, need to specify lambda (positive) and h (initial value)
  if (func == "decay") {
    return(sapply(X, decay <- function(x) {
      h * exp(-lambda * x)
    }))    
  }
  # if func is quadratic, need to specify coefficients a,b, and c
  if (func == "quadratic") {
    return(sapply(X, quadratic <- function(x) {
      a * x^2 + b* x + c
    }))
  }
  # you can try other funcs! 
}

## Obtain posterior distribution 
## Data given
D.train <- data.frame(x = c(1,2,3,4),
                     y = c(1,0.3,0.4,0.9))
D.test <- data.frame(x = c(1.9,2.3,7,10),
                     y = c(2,0.5,3,0.87))
D.test <- data.frame(x = 5,
                     y = 1.1)
# Step 1: Construct the kernel and mean
K.XtrainXtrain <- calcK(D.train$x, D.train$x, 
                        kernel = "SE", h = 0.3, lambda = 3)
K.XtrainXtest <- calcK(D.train$x, D.test$x, 
                       kernel = "SE", h = 0.3, lambda = 3)
K.XtestXtrain <- calcK(D.test$x, D.train$x, 
                       kernel = "SE", h = 0.3, lambda = 3)
K.XtestXtest <- calcK(D.test$x, D.test$x, 
                      kernel = "SE", h = 0.3, lambda = 3)
mean.train <- setMean(D.train$x)
mean.test <- setMean(D.test$x)

# Step 2: Get the posterior mean and cov matrix (eq. 3.8 and 3.9)
# do Cholesky decomposition for faster inversion
# <- chol(K.XtrainXtrain)
mean.pos <- mean.test + K.XtestXtrain %*% solve(K.XtrainXtrain) %*% (D.train$y - mean.train)
K.Xtest.pos <- K.XtestXtest - K.XtestXtrain %*% solve(K.XtrainXtrain) %*% t(K.XtestXtrain)

## Examples
# 1. Suppose we have a sequence of x's (times) with a mean 
#    function specified. We simulate the GP.
X <- seq(from = 0, to = 5, length = 1000)
K.XX <- calcK(X, X, "SE", h = 0.3, lambda = 0.3) 
mean.X <- setMean(X, func = "decay", lambda = 0.5, h = 10) 
# generate some sample paths 
n.samples <- 5
path_vals <- matrix(rep(0,length(X)*n.samples), 
                    ncol = n.samples)
for (i in 1:n.samples) {
  path_vals[,i] <- mvrnorm(1, mean.X, K.XX)
}
df <- cbind(x = X, 
            y = as.data.frame(path_vals))
df <- melt(df, id="x")
# compute the mean of sample paths 
mean.emp <- cbind(x = X, 
                  y = as.data.frame(rowMeans(path_vals)))
# plot graph of sample paths
ggplot(data = df, aes(x = x, y = value)) +
  geom_line(aes(group=variable, colour = variable),
            alpha = 0.8) +
  theme_classic() +
  geom_line(data = mean.emp, aes(x = x, y = `rowMeans(path_vals)`),
            lty = "dashed") +
  geom_line(data = as.data.frame(cbind(x = X, y = mean.X)),
            aes(x = x, y = `mean.X`),
            alpha = 0.8) +
  ylab("paths, f(x)") +
  theme(legend.position = "none") 

# 2. Suppose data is given and we want to estimate posterior mean   
X <- c(0,0.9,1.5,2,4)
y = c(0.1,-0.5,0.9,1,0.2)
K.XX <- calcK(X, X, "SE", h = 0.3, lambda = 0.3) 
mean.X <- setMean(X, func = "zero") 
# choose new X points
X.ast <- c(0.5,3,6,10)
# compute the covariance kernels 
K.XastXast <- calcK(X.ast, X.ast, "SE", h = 0.3, lambda = 0.3) 
K.XastX <- calcK(X.ast, X, kernel = "SE", h = 0.3, lambda = 0.3)
K.XXast <- calcK(X, X.ast, kernel = "SE", h = 0.3, lambda = 0.3)
mean.Xast <- setMean(X.ast, func = "zero")
# obtain the posterior mean and covariance
mean.pos <- mean.Xast + K.XastX %*% solve(K.XX) %*% (y - mean.X)
K.Xtest.pos <- K.XastXast - K.XastX %*% solve(K.XX) %*% t(K.XastX)

# compare graphs of sample paths with and w/o data
# first, do with two points observed
x <- c(seq(from = 0, to = 1.49, by = 0.01),
       seq(from = 1.51, to = 3.49, by = 0.01),
       seq(from = 3.51, to = 5, by = 0.01))
K.xx <- calcK(x, x, "SE", h = 0.3, lambda = 0.3)
mean.x <- setMean(x, func = "zero") 
f.obs <- data.frame(x = c(1.5,3.5),
           y = c(2,5)) 
mean.xobs <- setMean(f.obs$x, func = "zero") 
K.xobsxobs <- calcK(f.obs$x, f.obs$x, "SE", h = 0.3, lambda = 0.3)
K.xobsx <- calcK(f.obs$x, x, "SE", h = 0.3, lambda = 0.3)
K.xxobs <- calcK(x, f.obs$x, "SE", h = 0.3, lambda = 0.3)
# obtain posterior mean and cov kernel
mean.pos <- mean.x + K.xxobs %*% solve(K.xobsxobs) %*% (f.obs$y - mean.xobs)
K.x.pos <- K.xx - K.xxobs %*% solve(K.xobsxobs) %*% t(K.xxobs)

n.samples <- 50
path_vals <- matrix(rep(0,length(x)*n.samples), 
                    ncol = n.samples)
for (i in 1:n.samples) {
  path_vals[,i] <- mvrnorm(1, mean.pos, K.x.pos)
}
df <- cbind(x = x, 
            y = as.data.frame(path_vals))
#observed <- rbind(c(1.5,rep(2,times = 50)),
#                  c(3.5,rep(5,times = 50))) %>% as.data.frame()
#rbind(df, observed)
df <- melt(df, id="x")
# compute the mean of sample paths 
mean.emp <- cbind(x = x, 
                  y = as.data.frame(rowMeans(path_vals)))
plot2 <- ggplot(data = df, aes(x = x, y = value)) +
  geom_line(aes(group=variable, colour = variable),
            alpha = 0.8) +
  theme_classic() +
  geom_line(data = mean.emp, aes(x = x, y = `rowMeans(path_vals)`),
            lty = "dashed") +
  geom_line(data = as.data.frame(cbind(x = x, y = mean.pos)),
            aes(x = x, y = `mean.pos`),
            alpha = 0.8) +
  ylab(NULL) +
  theme(legend.position = "none") +
  geom_point(data = f.obs, aes(x = x, y = y),
             size = 2.5)

# next, do w/o data
x <- seq(from = 0, to = 5, by = 0.01)
K.xx <- calcK(x, x, "SE", h = 0.3, lambda = 0.3)
mean.x <- setMean(x, func = "zero") 
n.samples <- 50
path_vals <- matrix(rep(0,length(x)*n.samples), 
                    ncol = n.samples)
for (i in 1:n.samples) {
  path_vals[,i] <- mvrnorm(1, mean.x, K.xx)
}
df <- cbind(x = x, 
            y = as.data.frame(path_vals))
df <- melt(df, id="x")
# compute the mean of sample paths 
mean.emp <- cbind(x = x, 
                  y = as.data.frame(rowMeans(path_vals)))
plot1 <- ggplot(data = df, aes(x = x, y = value)) +
  geom_line(aes(group=variable, colour = variable),
            alpha = 0.8) +
  theme_classic() +
  geom_line(data = mean.emp, aes(x = x, y = `rowMeans(path_vals)`),
            lty = "dashed") +
  geom_line(data = as.data.frame(cbind(x = x, y = mean.x)),
            aes(x = x, y = `mean.x`),
            alpha = 0.8) +
  ylab("paths, f(x)") +
  theme(legend.position = "none") 
# now, plot!
gridExtra::grid.arrange(plot1,plot2,nrow = 1)