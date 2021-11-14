## Section First
# Preparation
N = 500
set.seed(1000) 
Z = runif(N, -2, 2) 
T = rbinom(N, 1,0.5)*2-1 
e = rnorm(N, 0,1/3)

Y = 2*Z*ifelse(T==1,1,0)+(Z-1)*ifelse(T==-1,1,0)+e

plot(Z*T/2,Y) 
const = rep(1, N) 
W = cbind(const, Z)

minimizer = matrix(nrow = 200, ncol = 200) 
for (i in 1:200){ 
  for(j in 1:200){ 
    gamma0 = c(i,j)/100 
    minimizer[i,j] = 1/500*sum((Y -T/2 * W %*% gamma0)ˆ2) 
  } 
  }

index = which(minimizer==min(minimizer),arr.ind=TRUE) 
print(index)

# model: Y ~ (1, Z)*T/2 
newW = W*T/2 
lm(Y ~ newW-1)#no intercept

## Section 2 Simulation
# %% [code]
# Numerical Studies

## SubSection 1 Continuous Responses

# For continuous responses, we generated $N$ independent Gaussian samples from a regression model.

# %% [code]
# introduce required package MASS and GLMnet
library(MASS)
library(glmnet)

# %% [code]
maineffect = function(p, basic = 6){
  parameter = c(1/sqrt(basic), 0, 0,rep(1/(2*sqrt(basic)),8),rep(0, p-10))
  return(parameter)
}

mat = function(p, rho){
  constant = rep(1, p)
  variance_covariance = (1 - rho)*diag(p) + rho*constant %*% t(constant)
  return(variance_covariance)
}

FullRegression = function(p,rho,basic =6){
  Cov = rep(0, 500)
  N = 100
  sigma0 = sqrt(2)
  alpha = 0.8
  gamma = c(0.4, 0.8, -0.8, 0.8, -0.8, rep(0, p+1-5))
  
  beta = maineffect(p, basic)
  
  variance_covariance = mat(p, rho)
  
  newN = 10000
  Z_generated = mvrnorm(newN, rep(0, p), variance_covariance)
  W_generated = cbind(rep(1, newN), Z_generated)
  delta  = 1.6*(0.5+Z_generated[,1]-Z_generated[,2]+Z_generated[,3]-Z_generated[,4]+Z_generated[,1]*Z_generated[,2])
  
  for(step in 1:500){
    e = rnorm(N)
    Z = mvrnorm(N, rep(0, p), variance_covariance)
    T = rbinom(N, 1, 0.5)*2-1
    const = rep(1, N)
    W = cbind(const,Z)
    Y = (W %*% beta)^2 + (W %*%gamma + alpha*Z[,1]*Z[,2])*T + e*sigma0
    covariate = cbind(W,W*T)
    
    dat = data.frame(cbind(Y, covariate))
    
    lambda_seq <- 10^seq(2, -2, by = -.1)
    
    lasso_model = cv.glmnet(covariate, dat$V1, alpha = 1, intercept = FALSE, nlambda = 40, 
                            family = "gaussian", standardize = FALSE)
    lambda_best = lasso_model$lambda.min
    lasso_best = glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, family = "gaussian", lambda = lambda_best, standardize = FALSE)
    coeff = coef(lasso_best)
    coeff = coeff[-1]
    coeff = coeff[-(1:(p+1))]
    
    if (!all(coeff==0)){
      score = W_generated %*% coeff
      if(length(unique(score))!=1){
        Cov[step] = cor(score, delta, method = "spearman")
      }
    }
  }
  return(Cov)
}

New = function(p,rho,basic =6){
  Cov = rep(0, 500)
  N = 100
  sigma0 = sqrt(2)
  alpha = 0.8
  gamma = c(0.4, 0.8, -0.8, 0.8, -0.8, rep(0, p+1-5))
  
  beta = maineffect(p, basic)
  
  variance_covariance = mat(p, rho)
  
  newN = 10000
  Z_generated = mvrnorm(newN, rep(0, p), variance_covariance)
  W_generated = cbind(rep(1, newN), Z_generated)
  delta  = 1.6*(0.5+Z_generated[,1]-Z_generated[,2]+Z_generated[,3]-Z_generated[,4]+Z_generated[,1]*Z_generated[,2])
  
  for(step in 1:500){
    e = rnorm(N)
    Z = mvrnorm(N, rep(0, p), variance_covariance)
    T = rbinom(N, 1, 0.5)*2-1
    const = rep(1, N)
    W = cbind(const,Z)
    Y = (W %*% beta)^2 + (W %*%gamma + alpha*Z[,1]*Z[,2])*T + e*sigma0
    
    
    Wstar = cbind(const, Z)*T/2
    
    dat = data.frame(cbind(Y, Wstar))
    lasso_model = cv.glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, nlambda = 20, family = "gaussian", standardize = FALSE)
    lambda_best = lasso_model$lambda.min
    lasso_best = glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, family = "gaussian", lambda = lambda_best, standardize = FALSE)
    coeff = coef(lasso_best)
    coeff = coeff[-1]
    
    
    if (!all(coeff==0)){
      score = W_generated %*% coeff
      if(length(unique(score))!=1){
        Cov[step] = cor(score, delta, method = "spearman")
      }
    }
  }
  return(Cov)
}

NewAugmented = function(p,rho,basic =6){
  Cov = rep(0, 500)
  N = 100
  sigma0 = sqrt(2)
  alpha = 0.8
  gamma = c(0.4, 0.8, -0.8, 0.8, -0.8, rep(0, p+1-5))
  
  beta = maineffect(p, basic)
  
  variance_covariance = mat(p, rho)
  
  newN = 10000
  Z_generated = mvrnorm(newN, rep(0, p), variance_covariance)
  W_generated = cbind(rep(1, newN), Z_generated)
  delta  = 1.6*(0.5+Z_generated[,1]-Z_generated[,2]+Z_generated[,3]-Z_generated[,4]+Z_generated[,1]*Z_generated[,2])
  
  for(step in 1:500){
    e = rnorm(N)
    Z = mvrnorm(N, rep(0, p), variance_covariance)
    T = rbinom(N, 1, 0.5)*2-1
    const = rep(1, N)
    W = cbind(const,Z)
    Y = (W %*% beta)^2 + (W %*%gamma + alpha*Z[,1]*Z[,2])*T + e*sigma0
    
    
    dat1 = data.frame(cbind(Y, W))
    lasso_model_xi = cv.glmnet(as.matrix(dat1[-1]), dat1$V1, alpha = 1, intercept = FALSE, nlambda = 20, family = "gaussian", standardize = FALSE)
    lambda_best = lasso_model_xi$lambda.min
    lasso_best = glmnet(as.matrix(dat1[-1]), dat1$V1, alpha = 1, intercept = FALSE, lambda = lambda_best, family = "gaussian", standardize = FALSE)
    xi = coef(lasso_best)
    xi = xi[-1]
    
    Wstar = cbind(const, Z)*T/2
    newY = Y -  W %*% xi
    
    dat = data.frame(cbind(newY, Wstar))
    lambda_seq <- 10^seq(2, -2, by = -.1)
    lasso_model = cv.glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, nlambda = 20, family = "gaussian", standardize = FALSE)
    lambda_best = lasso_model$lambda.min
    lasso_best = glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, family = "gaussian", lambda = lambda_best, standardize = FALSE)
    coeff = coef(lasso_best)
    coeff = coeff[-1]
    
    
    if (!all(coeff==0)){
      score = W_generated %*% coeff
      if(length(unique(score))!=1){
        Cov[step] = cor(score, delta, method = "spearman")
      }
    }
  }
  return(Cov)
}

# %% [code]
Cov = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov2 = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov3 = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov4 = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov5 = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov6 = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov7 = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov8 = cbind(rep(0, 500), rep(0, 500), rep(0, 500))

# %% [code]
Cov[,1] = FullRegression(50, 0, 6)
Cov[,2] = New(50, 0, 6)
Cov[,3] = NewAugmented(50, 0, 6)
Cov2[,1] = FullRegression(1000, 0, 6)
Cov2[,2] = New(1000, 0, 6)
Cov2[,3] = NewAugmented(1000, 0, 6)

# %% [code]
Cov3[,1] = FullRegression(50, 1/3, 6)
Cov3[,2] = New(50, 1/3, 6)
Cov3[,3] = NewAugmented(50, 1/3, 6)
Cov4[,1] = FullRegression(1000, 1/3, 6)
Cov4[,2] = New(1000, 1/3, 6)
Cov4[,3] = NewAugmented(1000, 1/3, 6)

# %% [code]
Cov5[,1] = FullRegression(50, 0, 3)
Cov5[,2] = New(50, 0, 3)
Cov5[,3] = NewAugmented(50, 0, 3)
Cov6[,1] = FullRegression(1000, 0, 3)
Cov6[,2] = New(1000, 0, 3)
Cov6[,3] = NewAugmented(1000, 0, 3)

# %% [code]
Cov7[,1] = FullRegression(50, 1/3, 3)
Cov7[,2] = New(50, 1/3, 3)
Cov7[,3] = NewAugmented(50, 1/3, 3)
Cov8[,1] = FullRegression(1000, 1/3, 3)
Cov8[,2] = New(1000, 1/3, 3)
Cov8[,3] = NewAugmented(1000, 1/3, 3)

# %% [code]
group1 = rep(1, 500)
group2 = rep(2, 500)
group3 = rep(3, 500)

situation2 = rep(2, 500)
situation1 = rep(1, 500)
x = c(group1, group1, group2, group2, group3, group3)
y = c(situation1, situation2, situation1, situation2, situation1, situation2)

# %% [code]
Covv = c(Cov2[,1],Cov[,1],Cov2[,2],Cov[,2], Cov2[,3], Cov[,3])
Covv2 = c(Cov4[,1],Cov3[,1],Cov4[,2],Cov3[,2], Cov4[,3], Cov3[,3])
Covv3 = c(Cov6[,1],Cov5[,1],Cov6[,2],Cov5[,2], Cov6[,3], Cov5[,3])
Covv4 = c(Cov8[,1],Cov7[,1],Cov8[,2],Cov7[,2], Cov8[,3], Cov7[,3])

# %% [code]
par(mfrow = c(2,2))
boxplot(Covv ~ y+x, dat = data.frame(Covv, x, y), names = c("Full", "", "New","","NewAu",""),xlab ="")
boxplot(Covv2  ~ y+x, dat = data.frame(Covv2, x, y), names = c("Full", "", "New","","NewAu",""),xlab ="")
boxplot(Covv3 ~ y+x, dat = data.frame(Covv3, x, y), names = c("Full", "", "New","","NewAu",""),xlab ="")
boxplot(Covv4 ~ y+x, dat = data.frame(Covv4, x, y), names = c("Full", "", "New","","NewAu",""),xlab ="")

# %% [code]
# SubSection 2 Binary Response

# %% [code]
maineffect = function(p, basic = 6){
  parameter = c(1/sqrt(basic), 0, 0,rep(1/(2*sqrt(basic)),8),rep(0, p-10))
  return(parameter)
}

mat = function(p, rho){
  rho = 0
  constant = rep(1, p)
  variance_covariance = (1 - rho)*diag(p) + rho*constant %*% t(constant)
  return(variance_covariance)
}

FullRegressionbin = function(p,rho,basic =6){
  Cov = rep(0, 500)
  N = 100
  sigma0 = sqrt(2)
  alpha = 0.8
  gamma = c(0.4, 0.8, -0.8, 0.8, -0.8, rep(0, p+1-5))
  
  beta = maineffect(p, basic)
  
  variance_covariance = mat(p, rho)
  
  newN = 10000
  Z_generated = mvrnorm(newN, rep(0, p), variance_covariance)
  W_generated = cbind(rep(1, newN), Z_generated)
  delta= pnorm(1/sigma0*((W_generated %*% beta)^2 + (W_generated %*%gamma + alpha*Z_generated[,1]*Z_generated[,2])))-pnorm(1/sigma0*((W_generated %*% beta)^2 - (W_generated %*%gamma + alpha*Z_generated[,1]*Z_generated[,2])))
  
  for(step in 1:500){
    e = rnorm(N)
    Z = mvrnorm(N, rep(0, p), variance_covariance)
    T = rbinom(N, 1, 0.5)*2-1
    const = rep(1, N)
    W = cbind(const,Z)
    Y = ifelse((W %*% beta)^2 + (W %*%gamma + alpha*Z[,1]*Z[,2])*T+e*sigma0>=0,1,0)
    covariate = cbind(W,W*T)
    
    dat = data.frame(cbind(Y, covariate))
    lambda_seq <- 10^seq(2, -2, by = -.1)
    lasso_model = cv.glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, nlambda = 20, 
                            family = "binomial")
    lambda_best = lasso_model$lambda.1se
    lasso_best = glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, family = "binomial", lambda = lambda_best, exact = T)
    coeff = coef(lasso_best)
    coeff = coeff[-1]
    coeff = coeff[-(1:(p+1))]
    
    if (!all(coeff==0)){
      score = W_generated %*% coeff
      if(length(unique(score))!=1){
        Cov[step] = cor(score, delta, method = "spearman")
      }
    }
  }
  return(Cov)
}

Newbin = function(p,rho,basic =6){
  Cov = rep(0, 500)
  N = 100
  sigma0 = sqrt(2)
  alpha = 0.8
  gamma = c(0.4, 0.8, -0.8, 0.8, -0.8, rep(0, p+1-5))
  
  beta = maineffect(p, basic)
  
  variance_covariance = mat(p, rho)
  
  newN = 10000
  Z_generated = mvrnorm(newN, rep(0, p), variance_covariance)
  W_generated = cbind(rep(1, newN), Z_generated)
  delta= pnorm(1/sigma0*((W_generated %*% beta)^2 + (W_generated %*%gamma + alpha*Z_generated[,1]*Z_generated[,2])))-pnorm(1/sigma0*((W_generated %*% beta)^2 - (W_generated %*%gamma + alpha*Z_generated[,1]*Z_generated[,2])))
  
  for(step in 1:500){
    e = rnorm(N)
    Z = mvrnorm(N, rep(0, p), variance_covariance)
    T = rbinom(N, 1, 0.5)*2-1
    const = rep(1, N)
    W = cbind(const,Z)
    Y = ifelse((W %*% beta)^2 + (W %*%gamma + alpha*Z[,1]*Z[,2])*T+e*sigma0>=0,1,0)
    covariate = cbind(W,W*T)
    
    Wstar = cbind(const, Z)*T/2
    
    dat = data.frame(cbind(Y, Wstar))
    lambda_seq <- 10^seq(2, -2, by = -.1)
    lasso_model = cv.glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, nlambda = 20, family = "binomial")
    lambda_best = lasso_model$lambda.min
    lasso_best = glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, family = "binomial", lambda = lambda_best)
    coeff = coef(lasso_best)
    coeff = coeff[-1]
    
    
    if (!all(coeff==0)){
      score = W_generated %*% coeff
      if(length(unique(score))!=1){
        Cov[step] = cor(score, delta, method = "spearman")
      }
    }
  }
  return(Cov)
}

NewAugmentedbin = function(p,rho,basic =6){
  Cov = rep(0, 500)
  N = 100
  sigma0 = sqrt(2)
  alpha = 0.8
  gamma = c(0.4, 0.8, -0.8, 0.8, -0.8, rep(0, p+1-5))
  
  beta = maineffect(p, basic)
  
  variance_covariance = mat(p, rho)
  
  newN = 10000
  Z_generated = mvrnorm(newN, rep(0, p), variance_covariance)
  W_generated = cbind(rep(1, newN), Z_generated)
  delta= pnorm(1/sigma0*((W_generated %*% beta)^2 + (W_generated %*%gamma + alpha*Z_generated[,1]*Z_generated[,2])))-pnorm(1/sigma0*((W_generated %*% beta)^2 - (W_generated %*%gamma + alpha*Z_generated[,1]*Z_generated[,2])))
  
  for(step in 1:500){
    e = rnorm(N)
    Z = mvrnorm(N, rep(0, p), variance_covariance)
    T = rbinom(N, 1, 0.5)*2-1
    const = rep(1, N)
    W = cbind(const,Z)
    Y = ifelse((W %*% beta)^2 + (W %*%gamma + alpha*Z[,1]*Z[,2])*T+e*sigma0>=0,1,0)
    
    
    dat1 = data.frame(cbind(Y, W))
    lambda_seq <- 10^seq(2, -2, by = -.1)
    lasso_model_xi = cv.glmnet(as.matrix(dat1[-1]), dat1$V1, alpha = 1, intercept = FALSE, lambda = lambda_seq, family = "binomial")
    xi = coef(lasso_model_xi, s = lasso_model_xi$lambda.min)
    xi = xi[-1]
    
    Wstar = cbind(const, Z)*T/2
    newY = Y -  W %*% xi
    
    dat = data.frame(cbind(newY, Wstar))
    #lambda_seq <- 10^seq(2, -2, by = -.1)
    lasso_model = cv.glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, nlambda = 20, family = "binomial")
    lambda_best = lasso_model$lambda.min
    lasso_best = glmnet(as.matrix(dat[-1]), dat$V1, alpha = 1, intercept = FALSE, family = "binomial", lambda = lambda_best)
    coeff = coef(lasso_best)
    coeff = coeff[-1]
    
    
    if (!all(coeff==0)){
      score = W_generated %*% coeff
      if(length(unique(score))!=1){
        Cov[step] = cor(score, delta, method = "spearman")
      }
    }
  }
  return(Cov)
}

# %% [code]
Covbin = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov2bin = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov3bin = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov4bin = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov5bin = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov6bin = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov7bin = cbind(rep(0, 500), rep(0, 500), rep(0, 500))
Cov8bin = cbind(rep(0, 500), rep(0, 500), rep(0, 500))

# %% [code]
Covbin[,1] = FullRegressionbin(50, 0, 6)
Covbin[,2] = Newbin(50, 0, 6)
Covbin[,3] = NewAugmentedbin(50, 0, 6)
Cov2bin[,1] = FullRegressionbin(1000, 0, 6)
Cov2bin[,2] = Newbin(1000, 0, 6)
Cov2bin[,3] = NewAugmentedbin(1000, 0, 6)

# %% [code]
Cov3bin[,1] = FullRegressionbin(50, 1/3, 6)
Cov3bin[,2] = Newbin(50, 1/3, 6)
Cov3bin[,3] = NewAugmentedbin(50, 1/3, 6)
Cov4bin[,1] = FullRegressionbin(1000, 1/3, 6)
Cov4bin[,2] = Newbin(1000, 1/3, 6)
Cov4bin[,3] = NewAugmentedbin(1000, 1/3, 6)

# %% [code]
Cov5bin[,1] = FullRegressionbin(50, 0, 3)
Cov5bin[,2] = Newbin(50, 0, 3)
Cov5bin[,3] = NewAugmentedbin(50, 0, 3)
Cov6bin[,1] = FullRegressionbin(1000, 0, 3)
Cov6bin[,2] = Newbin(1000, 0, 3)
Cov6bin[,3] = NewAugmentedbin(1000, 0, 3)

# %% [code]
Cov7bin[,1] = FullRegression(50, 1/3, 3)
Cov7bin[,2] = New(50, 1/3, 3)
Cov7bin[,3] = NewAugmented(50, 1/3, 3)
Cov8bin[,1] = FullRegression(1000, 1/3, 3)
Cov8bin[,2] = New(1000, 1/3, 3)
Cov8bin[,3] = NewAugmented(1000, 1/3, 3)

# %% [code]
group1 = rep(1, 500)
group2 = rep(2, 500)
group3 = rep(3, 500)

situation2 = rep(2, 500)
situation1 = rep(1, 500)
x = c(group1, group1, group2, group2, group3, group3)
y = c(situation1, situation2, situation1, situation2, situation1, situation2)

# %% [code]
Covvbin = c(Cov2bin[,1],Covbin[,1],Cov2bin[,2],Covbin[,2], Cov2bin[,3], Covbin[,3])
Covv2bin = c(Cov4bin[,1],Cov3bin[,1],Cov4bin[,2],Cov3bin[,2], Cov4bin[,3], Cov3bin[,3])
Covv3bin = c(Cov6bin[,1],Cov5bin[,1],Cov6bin[,2],Cov5bin[,2], Cov6bin[,3], Cov5bin[,3])
Covv4bin = c(Cov8bin[,1],Cov7bin[,1],Cov8bin[,2],Cov7bin[,2], Cov8bin[,3], Cov7bin[,3])

# %% [code]
par(mfrow = c(2,2))
boxplot(Covvbin ~ y+x, dat = data.frame(Covvbin, x, y), names = c("Full", "", "New","","NewAu",""),xlab ="")
boxplot(Covv2bin  ~ y+x, dat = data.frame(Covv2bin, x, y), names = c("Full", "", "New","","NewAu",""),xlab ="")
boxplot(Covv3bin ~ y+x, dat = data.frame(Covv3bin, x, y), names = c("Full", "", "New","","NewAu",""),xlab ="")
boxplot(Covv4bin ~ y+x, dat = data.frame(Covv4bin, x, y), names = c("Full", "", "New","","NewAu",""),xlab ="")
