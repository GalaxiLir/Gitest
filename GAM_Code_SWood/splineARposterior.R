# spline estimation under autoregressive error with a confidence interval/ a distribution 
# relative to the posterior distribution of the hyperparameter/smoothing parameter

# Install and load required packages
install.packages("mgcv")
install.packages("bayesm")
library(mgcv)
library(bayesm)

# Generate synthetic data with Autoregressive error structure
set.seed(123)
n <- 400
x <- seq(0, 1, length.out = n)
y <- sin(2 * pi * x)  + .7*arima.sim(model = list(ar = 0.7), n = n)
ystar<- sin(2 * pi * x)

# Fit a spline Generalized Additive Mixed Models model autoregressive error
model <- gamm(y ~ s(x), correlation = corAR1(form = ~1), data = data.frame(x = x, y = y))
model <- gamm(y ~ s(x), correlation = corAR1(), data = data.frame(x = x, y = y))
summary(model$gam)
summary(model$lme)
# Plot the results
plot(model$gam,pages=1)
lines(x,y-mean(y))
lines(x,ystar,col=2)
plot(model$lme)
anova(model$gam) 
gam.check(model$gam) # simple checking plots
## Raw residuals still show correlation, of course...
acf(residuals(model$gam),main="raw residual ACF")
## But standardized are now fine...
acf(residuals(model$lme,type="normalized"),main="standardized residual ACF")




# Bayesian estimation of the smoothing parameter
# Assuming a Gamma prior for the smoothing parameter
nIter <- 5000
burn <- 1000
thin <- 5
prior_shape <- 2
prior_rate <- 1

# Create a matrix to store the posterior samples
post_samples <- matrix(NA, nrow = nIter, ncol = 2)

# MCMC loop
for (i in 1:nIter) {
  model$sp <- rgamma(1, shape = prior_shape + 0.5 * length(model$residuals),
                     rate = prior_rate + 0.5 * sum(model$residuals^2))
  fit <- gamm(y ~ s(x), correlation = corAR1(form = ~1), data = data.frame(x = x, y = y),
              sp = model$sp)
  post_samples[i, ] <- c(fit$sp, sum(gam.check(fit)$residual.lk))
}

# Discard burn-in and thin the samples
post_samples <- post_samples[(burn + 1):(nIter - 1), seq(1, ncol(post_samples), by = thin)]

# Plot posterior distribution of the smoothing parameter
hist(post_samples[, 1], main = "Posterior Distribution of Smoothing Parameter",
     xlab = "Smoothing Parameter", col = "lightblue", border = "black")

# Confidence interval for the smoothing parameter
quantiles <- quantile(post_samples[, 1], c(0.025, 0.975))
cat("95% Confidence Interval for Smoothing Parameter:", quantiles, "\n")