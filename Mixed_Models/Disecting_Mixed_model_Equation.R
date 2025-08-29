## R Script: Disecting Mixed Model Equation
## Set working Directory

  setwd("~/Documents/Mixed_Models")

## ----message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------------
library(agridat) 
  library(lme4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Load the Australia.soybean data from agridat package
    australia.soybean<-read.csv(file="./Data/australia.soybean.csv")

# Get the structure of data
    str(australia.soybean)
# Convert env, loc and year into factors
    australia.soybean$env<-as.factor(australia.soybean$env)
    australia.soybean$loc<-as.factor(australia.soybean$loc)
    australia.soybean$year<-as.factor(australia.soybean$year)
    australia.soybean$gen<-as.factor(australia.soybean$gen)
# Data can also be upload directly from agridat package
    data(australia.soybean, package = "agridat")
    head(australia.soybean)
    str(australia.soybean)
# Subset the data (environment, genotypes and yield) now for mixed model equation
    demo.data<-australia.soybean[, c(1,4,5)]
# Look for number of environments.
    table(demo.data$env)
# Look for number of environments.
    table(demo.data$gen)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Building X matrix for fixed component, i.e., for environments
    X <- model.matrix(~+env, demo.data) # +env means I am adding In 
    head(X)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Building Z matrix for random component, i.e., for genotypes
    Z <- model.matrix(~-1 + gen, demo.data) # -1 here is for removing intercept


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Yield variable as response variable
    y <-demo.data$yield


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Cross product and transpose
    XtX <- t(X) %*% X
    XtX


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Cross product and transpose
    XtZ <- t(X) %*% Z
    XtZ[, 1:10]


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Cross product and transpose
    ZtX <- t(Z) %*% X
    ZtX[1:8, ]



## --------------------------------------------------------------------------------------------------------------------------------------------------
# Cross product and transpose
    ZtZ <- t(Z) %*% Z
    ZtZ[1:7,1:7]


## --------------------------------------------------------------------------------------------------------------------------------------------------
    sigmau <- 0.199 # Genotype variance
    sigmae <- 0.25 # Residual variance
    I <- diag(ncol(Z))  # assuming G = I, No fancy markers
    lambda <- sigmae/sigmau # shrinkage or Heritability


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Cross product and transpose
    Xty <- t(X) %*% y
    Xty


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Cross product and transpose
    Zty <- t(Z) %*% y
    head(Zty)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Left hand side 
    LHS1 <- cbind(XtX, XtZ)
    LHS2 <- cbind(ZtX, ZtZ + I * lambda)
    LHS <- rbind(LHS1, LHS2)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Right Hand side
    RHS <- rbind(Xty, Zty)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Sol function to get solution
    sol <- solve(LHS, RHS)
    dim(sol)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# First eight are fixed effects (BLUEs) for environments
  Blues.env<-data.frame(Blues.env=sol[1:8, ])
  Blues.env
  str(Blues.env)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Nine to rest are eight are random effects (BLUPs)
# BLUP
  Blups.gy<-sol[9:66, ]
  head(Blups.gy)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Final genotypic values/breeding values for grain yield
    blup.final<-data.frame(Yield.blups=Blups.gy+sol[1,1])


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Fit mixed model
    mixemodel.fit<-lmer(yield ~ env + (1 | gen), data = demo.data)
    summary(mixemodel.fit)
# Now extract the fixed effects (BLUEs) for environmenet
    Blues.env.lme4<-data.frame(Blues.env=mixemodel.fit@beta)
    Blues.env.lme4
# Now extract the random effects (BLUPs)
    Blups.gy.lme4<-ranef(mixemodel.fit)$gen
# Now extract the intercept and add it to random effects
    intercept <- fixef(mixemodel.fit)[1]
# Now add intercept to blups to get genotypic values
    Final_blups.gy.lme4<-data.frame(Yield.blup=intercept+ranef(mixemodel.fit)$gen)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Blues from both
    table(round(Blues.env$Blues.env,2)==round(Blues.env.lme4$Blues.env,2))
# Now for Blups
    table(round(blup.final$Yield.blups, 2)==round( Final_blups.gy.lme4$X.Intercept., 2))

