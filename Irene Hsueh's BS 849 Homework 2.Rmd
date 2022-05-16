---
title: "Irene Hsueh's BS 849 Homework 2"
author: "Irene Hsueh"
date: "2/7/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(rjags)
library(coda)
library(formatR)
```


# Frequentist Crude Logistic Regression 
```{r}
crude_frequency <- array(c(441, 423, 58, 147),
                         dim = c(2, 2), 
                         dimnames = list(disease = c("No MI", "MI"),
                                         exposure = c("<5 Cups", "5+ Cups")))
contingency_table_crude <- as.data.frame(as.table(crude_frequency))
crude_logistic_model1 <- glm(disease ~ exposure, weights=Freq, 
                             data=contingency_table_crude, family="binomial")
summary(crude_logistic_model1)


exposure_crude <- c(rep(1, 147), rep(0, 423), rep(1, 58), rep(0, 441))
disease_crude <- c(rep(1, 570), rep(0, 499))
table(disease_crude, exposure_crude)
crude_logistic_model2 <- glm(disease_crude ~ exposure_crude, family="binomial")
summary(crude_logistic_model2)
```



# Bayesian Logistic Regression 
```{r}
crude_model_bugs <- 
"model {

#Data Model
for(i in 1:N){
disease_crude[i] ~ dbin(p[i], 1)
logit(p[i]) <- beta0 + beta1*exposure_crude[i]
}
OR <- exp(beta1)

#Prior Distribution 
beta0 ~ dnorm(0, 0.0001)
beta1 ~ dnorm(0, 0.0001)
}"

crude_data <- list(N=1069, exposure_crude=exposure_crude, disease_crude=disease_crude)
crude_odds_model <- jags.model(textConnection(crude_model_bugs), data=crude_data, n.adapt=2000)
crude_model_gibbs <- update(crude_odds_model, n.iter=5000)

crude_model_test <- coda.samples(crude_odds_model, c("OR", "beta0", "beta1"), n.iter=10000)
summary(crude_model_test)
autocorr.plot(crude_model_test)
```



# Frequentist Adjusted Logistic Regression 
```{r}
adjusted_frequency <- array(c(327, 207, 33, 30, 114, 216, 25, 117),
                            dim = c(2, 2, 2), 
                            dimnames = list(disease = c("No MI", "MI"),
                                            exposure = c("<5 Cups", "5+ Cups"),
                                            smoking = c("Nonsmoker", "Smoker")))
adjusted_contingency_table <- as.data.frame(as.table(adjusted_frequency))
adjusted_logistic_model1 <- glm(disease ~ exposure + smoking, weights=Freq,
                                data=adjusted_contingency_table, family="binomial")
summary(adjusted_logistic_model1)


exposure_adjusted <- rep(c(1,0,1,0,1,0,1,0), c(30, 207, 117, 216, 33, 327, 25, 114)) 
smoking_adjusted <- rep(c(0, 1, 0, 1), c(237, 333, 360, 139))
disease_adjusted <- rep(c(1, 0), c(570, 499))
adjusted_logistic_model2 <- glm(disease_adjusted ~ exposure_adjusted + smoking_adjusted,
                                family="binomial")
summary(adjusted_logistic_model2)
```



# Bayesian Adjusted Logistic Regression 
```{r}
#step function is 1 if quantity ≥ 0, 0 if quantity < 0
#null_odds is 1 if OR>1, 0 if OR≤1 

adjusted_model_bugs <- 
"model {

#Data Model
for(i in 1:N){
disease_adjusted[i] ~ dbin(p[i], 1)
logit(p[i]) <- beta0 + beta1*exposure_adjusted[i] + beta2*smoking_adjusted[i]
}
OR <- exp(beta1)
null_odds <- step(OR-1)

#Prior Distribution 
beta0 ~ dnorm(0, 0.0001)
beta1 ~ dnorm(0, 0.0001)
beta2 ~ dnorm(0, 0.0001)
}"

adjusted_data <- list(N=1069, exposure_adjusted=exposure_adjusted,
                      disease_adjusted=disease_adjusted, smoking_adjusted=smoking_adjusted)
adjusted_odds_model <- jags.model(textConnection(adjusted_model_bugs), 
                                  data=adjusted_data, n.adapt=1000)
adjusted_model_gibbs <- update(adjusted_odds_model, n.iter=5000)

adjusted_model_test <- coda.samples(adjusted_odds_model, c("OR", "null_odds","beta0", "beta1", "beta2"), n.iter=10000)
summary(adjusted_model_test)
autocorr.plot(adjusted_model_test)
```



# Bayesian Stratified Logistic Regression for Interaction
```{r}
#step function is 1 if quantity ≥ 0, 0 if quantity < 0
#null_interaction is 1 if beta3>0, 0 if beta3≤0 

stratified_model_bugs <- 
"model {

#Data Model
for(i in 1:N){
disease_adjusted[i] ~ dbin(p[i], 1)
logit(p[i]) <- beta0 + beta1*exposure_adjusted[i] + beta2*smoking_adjusted[i] + beta3*exposure_adjusted[i]*smoking_adjusted[i]
}
OR_nosmoking <- exp(beta1)
OR_smoking <- exp(beta1 + beta3)
ROR <- exp(beta3)
null_interaction <- step(beta3-0)

#Prior Distribution 
beta0 ~ dnorm(0, 0.0001)
beta1 ~ dnorm(0, 0.0001)
beta2 ~ dnorm(0, 0.0001)
beta3 ~ dnorm(0, 0.0001)
}"

stratified_data <- list(N=1069, exposure_adjusted=exposure_adjusted,
                        disease_adjusted=disease_adjusted, smoking_adjusted=smoking_adjusted)
stratified_odds_model <- jags.model(textConnection(stratified_model_bugs),
                                    data=stratified_data, n.adapt=1000)
confounder_model_gibbs <- update(stratified_odds_model, n.iter=5000)

stratified_model_test <- coda.samples(stratified_odds_model, c("OR_nosmoking", "OR_smoking", "ROR", "null_interaction", "beta0", "beta1", "beta2", "beta3"), n.iter=10000)
summary(stratified_model_test)
autocorr.plot(stratified_model_test)
```



