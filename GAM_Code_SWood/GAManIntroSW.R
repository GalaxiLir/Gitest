#Generalised additive model an introduction S Wood
#debut au chap 2 LMM
require(gamair)


data(stomata)

m1 <- lm(area ~ CO2 + tree, stomata)
summary(m1)

m0 <- lm(area ~ CO2, stomata)
summary(m0)
anova(m0,m1)

m2 <- lm(area ~ tree, stomata)
anova(m2,m1)
#conclusion on ne peut pas conclure à un effet du CO2 ou des arbres qui soit significatif 
#etant donné que les données sont dépendantes de l'arbre en question
#seul le model contenant les deux effet donne de bon resultat selon fisher 
#model mixte



st <- aggregate(data.matrix(stomata), by=list(tree=stomata$tree),mean)#Matrice contenant la moyenne pour chaque arbre
st$CO2 <- as.factor(st$CO2);st
m3 <- lm(area ~ CO2, st)

summary(m3)
summary(m3)$sigma^2 - summary(m1)$sigma^2/4


library(nlme) # load nlme ‘library’, which contains data
data(Rail) # load data
Rail
m1 <- lm(travel ~ Rail,Rail)
 anova(m1)
 rt <- aggregate(data.matrix(Rail),by=list(Rail$Rail),mean) # average over Rail effect
 rt
 
 m0 <- lm(travel ~ 1,rt) # fit model to aggregated data
 sig <- summary(m1)$sigma # sig^2 is resid. var. component
  sigb <- (summary(m0)$sigma^2 - sig^2/3)^0.5
  # sigb^2 is the variance component for rail
  sigb
  sig
  
  
  library(nlme)
  data(Machines)
  names(Machines)

  attach(Machines) # make data available without ‘Machines$’
  interaction.plot(Machine,Worker,score)
  m1 <- lm(score ~ Worker*Machine,Machines)
  m0 <- lm(score ~ Worker + Machine,Machines)
  anova(m0,m1)
  summary(m1)$sigma^2
  Mach <- aggregate(data.matrix(Machines),by=
                      list(Machines$Worker,Machines$Machine),mean)
  Mach$Worker <- as.factor(Mach$Worker)
  Mach$Machine <- as.factor(Mach$Machine)
  m0 <- lm(score ~ Worker + Machine,Mach)
  anova(m0) 
  
 # LMM In general
# hypothèse  
  
  
  llm <- function(theta,X,Z,y) {
    ## untransform parameters...
    sigma.b <- exp(theta[1])
    sigma <- exp(theta[2])
    ## extract dimensions...
    n <- length(y); pr <- ncol(Z); pf <- ncol(X)
    ## obtain \hat \beta, \hat b...
    X1 <- cbind(X,Z)
    ipsi <- c(rep(0,pf),rep(1/sigma.b^2,pr))
    b1 <- solve(crossprod(X1)/sigma^2+diag(ipsi),
                t(X1)%*%y/sigma^2)
    ## compute log|Z’Z/sigma^2 + I/sigma.b^2|...
    ldet <- sum(log(diag(chol(crossprod(Z)/sigma^2 +
                                diag(ipsi[-(1:pf)])))))
    ## compute log profile likelihood...
    l <- (-sum((y-X1%*%b1)^2)/sigma^2 - sum(b1^2*ipsi) -
            n*log(sigma^2) - pr*log(sigma.b^2) - 2*ldet - n*log(2*pi))/2
    attr(l,"b") <- as.numeric(b1) ## return \hat beta and \hat b
    -l
  }
  
  
 