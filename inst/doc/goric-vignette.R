## ---- vinylidene, warning=FALSE, message=FALSE---------------------------
library(goric)
data(vinylidene)

## ---- lm-----------------------------------------------------------------
X <- model.matrix(~ dose, data=vinylidene)[,-1]
m <- lm(SDH ~ X, data=vinylidene)
coefficients(m)

## ---- mselset------------------------------------------------------------
m0 <- lm(SDH ~ 1, data=vinylidene)
m1 <- lm(SDH ~ X[,-1], data=vinylidene)
m2 <- lm(SDH ~ X[,-2], data=vinylidene)
m3 <- lm(SDH ~ X[,-3], data=vinylidene)
m4 <- lm(SDH ~ X[,-c(1,2)], data=vinylidene)
m5 <- lm(SDH ~ X[,-c(2,3)], data=vinylidene)
m6 <- lm(SDH ~ X[,-c(1,3)], data=vinylidene)

## ---- AICsel-------------------------------------------------------------
AIC(m, m0,m1,m2,m3,m4,m5,m6)

## ---- m4coef-------------------------------------------------------------
coefficients(m4)

