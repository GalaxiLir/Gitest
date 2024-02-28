  
##2D ploting function on agrid 
#x,y are the are the represent the discretisationn used on x,y from which the grid has been made
#Y is the vector containing the value of the function
plot_function_contour <- function(x,y, Y) {
  # Reshape Y to match the dimensions of X
  Z <- matrix(Y, nrow = length(unique(x)), ncol = length(unique(y)), byrow = FALSE)
  
  # Create a 2D contour plot
  contour(x, y, Z, main = "Function Evaluated on X", xlab = "X", ylab = "Y", nlevels = 20, col = rainbow(20))
   
  # Optionally, you can add a color scale legend
  color.legend(x = "topright", legend = "Y", fill = rainbow(40), cex = 0.8)
  
  # You can customize the plot further as per your preference
}
### Creating a data set to make inference on wether an additive multiplicative or
### full modelisation is requiered 
#testing if i can push
x=seq(f=0,t=1,l=100)
y=seq(f=0,t=1,l=100)

#takes all the possible products between x and y
X <- expand.grid(x = x, y = y)


addFun=function(x=NULL,y=NULL){
  X <- expand.grid(x = x, y = y)
  cos(2*2*pi*X[,1])+cos(2*2*pi*X[,2])
  #cos(2*2*pi*(X[,1]+4*X[,2]))
}

addFun <- function(x = NULL, y = NULL, fun1, fun2) {
  X <- expand.grid(x = x, y = y)
  Z1 <- fun1( X[, 1],  X[, 2])
  Z2 <- fun2( X[, 1],  X[, 2])
  Z1 + Z2
}

multFun <- function(x = NULL, y = NULL,fun0, fun1, fun2) {
  X <- expand.grid(x = x, y = y)
  Z0 <- fun0( X[, 1],  X[, 2])
  Z1 <- fun1( X[, 1],  X[, 2])
  Z2 <- fun2( X[, 1],  X[, 2])
  Z0+Z1*Z2
}

nonSepFun <- function(x = NULL, y = NULL, fun3) {
  X <- expand.grid(x = x, y = y)
  Z <- fun3( X[, 1],  X[, 2])
}
# Example functions
fun0 <- function(x, y) cos(2 * pi *x)
fun1 <- function(x, y) y*(y-1/8)*(y-3/8)*(y-1/4)*3+0.3*y
fun2 <- function(x, y) sin(2 * pi *x)
fun3 <- function(x, y) cos(2*pi*(2*x+y))
#Plot Fun1
plot(y,fun1(1,y),type = "l")
#plot Fun2 Fun0
plot(c(x,x),c(fun0(x,1),fun2(x,1)),type="n")
lines(x,fun0(x,1),col=2)
lines(x,fun2(x,1),col=4)



Yadd=addFun(x,y,fun1,fun2)
Ymul=multFun(x,y,fun0,fun1,fun2)
YnonSep=nonSepFun(x,y,fun3)
Yfull=Ymul+YnonSep


plot_function_contour(x,y,Yadd)
plot_function_contour(x,y,Ymul)
plot_function_contour(x,y,YnonSep)
plot_function_contour(x,y,Yfull)

### Error generation

n=length(Yadd)
#white noise
eps=rnorm(n=n)

# Autoregressive parameter
phi <- 0.8


# Generate autocorrelated errors using an AR(1) model
autocorrelated_errors <- arima.sim(model = list(ar = phi), n = n, innov = eps)

# Plot autocorrelated errors
plot(autocorrelated_errors, type = 'l', main = "Autocorrelated Errors", ylab = "Error")

acf(autocorrelated_errors)
pacf(autocorrelated_errors)

### Gam
library(mgcv)

## Additive data structure

#changing the signal to noise ratio manually
signal2noiseRinit=var(Yadd)/var(eps)#initial signal to noise ratio 
signal2noiseR=9 #target for signal to noise ratio 
sigma=sqrt(signal2noiseRinit/signal2noiseR)
var(Yadd)/var(sigma*eps) #verification

#main data.frame
Dat_F=as.data.frame(cbind(Yadd+sigma*eps,X))
colnames(Dat_F)=c("Y","x","y" )



### Additive framework


# Fit GAM model
model_additive <- gam(Y ~ s(x, bs = "cc",k = 25) + s(y,k = 25),data = Dat_F, method = "REML")

# Summary of model
summary(model_additive)

### Multiplicative Term Model

## Deriving a value for g
#from preceding additive model
g_y_pred <- predict(model_additive, type = "terms", terms = "s(y)")
# #or Fit a new GAM model for g(y)
# model_g_y <- gam(Y ~ s(y),data = Dat_F, method = "REML")
# # Predict g(y) values
# g_y_pred <- predict(model_g_y, type = "response")-model_g_y$coefficients[1]

# Fit GAM model with multiplicative term
model_multiplicative <- gam(Y ~  s(x) + g_y_pred+s(x, by = g_y_pred,bs="cc"),data = Dat_F, method = "REML")
multi_pred=predict(model_multiplicative, type = "response")
# Summary of model
summary(model_multiplicative)
plot(model_multiplicative)
### No constaints tensor product Splines Model Y=phi(x,y)
Dat_F=cbind(Dat_F,multi_pred)
# Fit GAM model with periodic constraint on x
model_te_splines <- gam(Y ~ multi_pred+te(x, y, bs = c("cc", "cr")),data = Dat_F, method = "REML")

plot(model_te_splines )
# Summary of model
summary(model_te_splines)

#testing hypothesis 

anova(model_additive, model_te_splines, test = "Chisq")
anova(model_additive, model_g_y, test = "Chisq")

anova(model_additive, model_multiplicative, test = "Chisq")
anova.gam(model_additive, model_multiplicative, test = "Chisq")
lrtest(model_additive, model_multiplicative)

#computing GLRT manually
dev_model_additive <- deviance(model_additive)
dev_model_multiplicative <- deviance(model_multiplicative)
dev_diff <- dev_model_additive - dev_model_multiplicative

p_value <- pchisq(dev_diff, df = 1, lower.tail = FALSE)
print(p_value)

### test to see if inputing directly s(y) directly changes deviance
# for signal2noiseR=9 différence in deviance 0.00165
model_dev1 <- gam(Y ~  s(x, bs = "cc",k = 25) + s(y,k = 25),data = Dat_F)
g_y_pred <- predict(model_dev1, type = "terms", terms = "s(y)")
model_dev2 <- gam(Y ~  s(x, bs = "cc",k = 25) + g_y_pred,data = Dat_F)
deviance(model_dev1)
deviance(model_dev2)

sum((residuals(model_dev1))^2)
sum((residuals(model_dev2))^2)
plot(model_dev2$residuals-model_dev1$residuals)
abline(h=0)
summary(model_dev2)
