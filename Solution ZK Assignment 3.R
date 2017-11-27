#### Solution ZK Assignment 3 ####
# clear environment
remove(list = ls())
library(psych)
library(lm.beta)
library(lattice)
library(r2glmm)
library(lmtest)
library(car)
library(cAIC4)
library(influence.ME)
library(lattice)
library(reshape2)
library(ggplot2)
library(HLMdiag)

#read in data
home_sample_3 = read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class/master/home_sample_3.csv") 

# copy data
home_sample_31 = home_sample_3
View(home_sample_31)

# Check the structure of that data frame
str(home_sample_3)

# asign ID as factors
home_sample_3$ID = factor(home_sample_3$ID)

# histograms of pain over time
hist(home_sample_3$pain1, breaks =10)
hist(home_sample_3$pain2, breaks =10)
hist(home_sample_3$pain3, breaks =10)
hist(home_sample_3$pain4, breaks =10)

# correlation of repeated variables
repeated_var = c("pain1", "pain2", "pain3", "pain4")
cor(home_sample_3[,repeated_var])


# convert data to long format
library(reshape2)

home_sample_3_long <- melt(home_sample_31, id.vars=c("ID","sex","age","STAI_trait","pain_cat","cortisol_serum", "cortisol_saliva", "mindfulness", "weight"), variable.name= "pain_time", value.name="pain_value")
View(home_sample_3_long)

# Rename and retype the time variable based on equal to day of pain report
names(home_sample_3_long)
home_sample_3_long$pain_time <- as.character(home_sample_3_long$pain_time)
mode(home_sample_3_long$pain_time)
home_sample_3_long$pain_time[home_sample_3_long$pain_time=="pain1"] <- 1
home_sample_3_long$pain_time[home_sample_3_long$pain_time=="pain2"] <- 2
home_sample_3_long$pain_time[home_sample_3_long$pain_time=="pain3"] <- 3
home_sample_3_long$pain_time[home_sample_3_long$pain_time=="pain4"] <- 4
home_sample_3_long$pain_time <- as.numeric(home_sample_3_long$pain_time)
mode(home_sample_3_long$pain_time)

# histograms
hist(home_sample_3_long$age, breaks=20)
hist(home_sample_3_long$weight, breaks=20)
hist(home_sample_3_long$STAI_trait, breaks=20)
hist(home_sample_3_long$pain_cat, breaks=20)
hist(home_sample_3_long$mindfulness, breaks=20)
hist(home_sample_3_long$cortisol_serum, breaks=20)

# descriptive data
describe(home_sample_3_long)
#  --- look 

# general linear model  
lm_fixed=lm(pain_value ~ age+sex+weight+STAI_trait+pain_cat+mindfulness+cortisol_serum+pain_time, data = home_sample_3_long)
summary(lm_fixed)  

# linear mixed model with random intercept for ID
library(lme4)
lme_1 <-lmer(pain_value~1+age+sex+weight+STAI_trait+pain_cat+mindfulness+cortisol_serum+pain_time+(1|ID),data=home_sample_3_long)
summary(lme_1)

# linear mixed model with random intercept AND slope for ID
lme_2 <-lmer(pain_value~1+age+sex+weight+STAI_trait+pain_cat+mindfulness+cortisol_serum+pain_time+(1+pain_time|ID),data=home_sample_3_long)
summary(lme_2)                           

# plot pain over time for each subject
# 1. plot with normal regression line
xyplot(cbind(home_sample_3_long$pain_value, home_sample_3_long$pred.rt.lme1)~pain_time|ID, xlab="time", ylab="pain value", main= "Predicted values connected with regression line", layout=c(4,5), type = c("p","r"), data=home_sample_3_long)

# 2. plot with regresssion line of random intercept model
xyplot(pain_value + fitted(lme_1)~pain_time|ID, data=home_sample_3_long,
      type= c("p", "l"), distribute.type = TRUE,
       xlab="time", ylab="pain value", main= "Random intercept regression line",
       layout=c(4,5))

# 3. plot with regression line of random intercept and random slope model
xyplot(pain_value + fitted(lme_2)~pain_time|ID, data=home_sample_3_long,
       type= c("p", "l"), distribute.type = TRUE,
       xlab="time", ylab="pain value", main= "Random intercept and random slope regression line",
       layout=c(4,5))

# model comparing, cAIC as Steinian type estimator based on adapted code from: Greven, S. and Kneib T. (2010) On the behaviour of marginal and conditional AIC in linear mixed models. Biometrika 97(4), 773-789.
library("cAIC4")
cAIC(lme_1)
cAIC(lme_2)
summary(lme_1)
summary(lme_2)
anova(lme_1, lme_2)

# Calculate marginal R squared by Nakagawa and Schielzeth (extended by Johnson), Article: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4368045/
r2beta(lme_1, method = "nsj") # if the CI doesnt overlap with zero, its good
r2beta(lme_2, method = "nsj")

# new model with quadratic time
lme_3 <-lmer(pain_value~1+age+sex+weight+STAI_trait+pain_cat+mindfulness+cortisol_serum+pain_time+I(pain_time^2)+(1+pain_time|ID),data=home_sample_3_long)
summary(lme_3)

# model characteristics
cAIC(lme_3)

# Calculate marginal R squared by Nakagawa and Schielzeth (extended by Johnson), Article: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4368045/
r2beta(lme_3, method = "nsj")
cAIC(lme_3)

# plot lme3 with regression line
home_sample_3_long$pred.pain.lme3<-predict(lme_3)
lme_3_plot <- xyplot(pain_value + fitted(lme_3)~pain_time|ID, data=home_sample_3_long,
                     type= c("p", "l"), distribute.type = TRUE,
                     xlab="time", ylab="pain value", main= "Model with quadratic term of time",
                     layout=c(4,5))
lme_3_plot

# save lme3 jpg
jpeg(filename="~/Dropbox/HU Berlin/3. Semester Lund University/PSYP13/Home Assignments/ex3_lme3_plot.jpg")
plot(lme_3_plot)
dev.off()

# function to extract standardized beta coefficients from mer models, such as one produced by lmer
# source: https://stackoverflow.com/questions/25142901/standardized-coefficients-for-lmer-model

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

stdCoef.merMod(lme_3)
confint(lme_3) # error messages are ok here, # no overelap with zero -> assumption ok!

# checking model-assumptions of lme3
# checking for influential outliers
influence(lme_3, group = "ID")$alt.fixed
influence(lme_3, obs = T)$alt.fixed

ggplot(data.frame(lev=hatvalues(lme_3),pearson=residuals(lme_3,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()


# linearity assumption
pred_lme_3 <- predict(object = lme_3) # predicted values against actual values
plot(x = pred_lme_3, y = home_sample_3_long$pain_value, 
      xlab = "Fitted Values", ylab = "Observed Values")
plot(lme_3, which = 1) # predicted values against residuals

# linearity of each predictor
predictors=c("age")
for(i in 1:length(predictors)){
  predictor_to_test = home_sample_3_long[,predictors[i]]
  print(
    ggplot(data.frame(x = predictor_to_test,pearson=residuals(lme_3,type="pearson")),
           aes(x=x,y=pearson)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      theme_bw() )
}


predictors=c("weight")
for(i in 1:length(predictors)){
  predictor_to_test = home_sample_3_long[,predictors[i]]
  print(
    ggplot(data.frame(x = predictor_to_test,pearson=residuals(lme_3,type="pearson")),
           aes(x=x,y=pearson)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      theme_bw() )
}

predictors=c("STAI_trait")
for(i in 1:length(predictors)){
  predictor_to_test = home_sample_3_long[,predictors[i]]
  print(
    ggplot(data.frame(x = predictor_to_test,pearson=residuals(lme_3,type="pearson")),
           aes(x=x,y=pearson)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      theme_bw() )
}

predictors=c("pain_cat")
for(i in 1:length(predictors)){
  predictor_to_test = home_sample_3_long[,predictors[i]]
  print(
    ggplot(data.frame(x = predictor_to_test,pearson=residuals(lme_3,type="pearson")),
           aes(x=x,y=pearson)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      theme_bw() )
}

predictors=c("mindfulness")
for(i in 1:length(predictors)){
  predictor_to_test = home_sample_3_long[,predictors[i]]
  print(
    ggplot(data.frame(x = predictor_to_test,pearson=residuals(lme_3,type="pearson")),
           aes(x=x,y=pearson)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      theme_bw() )
}

predictors=c("cortisol_serum")
for(i in 1:length(predictors)){
  predictor_to_test = home_sample_3_long[,predictors[i]]
  print(
    ggplot(data.frame(x = predictor_to_test,pearson=residuals(lme_3,type="pearson")),
           aes(x=x,y=pearson)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      theme_bw() )
}


# homoscedasticty assumption (homogeneity of variance)
plot(lme_3)
summary(lm(residuals(lme_3)^2 ~ home_sample_3_long[,"ID"])) # Levens model for heteroscedasticity -> not sig. 
# alternative way to compare 
home_sample_3_long$lme3_res<- residuals(lme_3) #extracts the residuals and places them in a new column in our original data table
home_sample_3_long$lme3_res_abs <-abs(home_sample_3_long$lme3_res) #creates a new column with the absolute value of the residuals
home_sample_3_long$lme3_res_abs_2 <- home_sample_3_long$lme3_res_abs^2 #squares the absolute values of the residuals to provide the more robust estimate
Levene.Model.F <- lm(lme3_res_abs_2 ~ ID, data=home_sample_3_long) #ANOVA of the squared residuals
anova(Levene.Model.F) #displays the results -> not sig.

# test multicollinearity
pairs.panels(home_sample_3_long, col = "red", lm = T)

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

kappa.mer <- function (fit,
                       scale = TRUE, center = FALSE,
                       add.intercept = TRUE,
                       exact = FALSE) {
  X <- fit@pp$X
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  nrp <- sum(1 * (nam == "(Intercept)"))
  if (nrp > 0) {
    X <- X[, -(1:nrp), drop = FALSE]
    nam <- nam[-(1:nrp)]
  }
  
  if (add.intercept) {
    X <- cbind(rep(1), scale(X, scale = scale, center = center))
    kappa(X, exact = exact)
  } else {
    kappa(scale(X, scale = scale, center = scale), exact = exact)
  }
}

colldiag.mer <- function (fit,
                          scale = TRUE, center = FALSE,
                          add.intercept = TRUE) {
  ## adapted from perturb::colldiag, method in Belsley, Kuh, and
  ## Welsch (1980).  look for a high condition index (> 30) with
  ## more than one high variance propotion.  see ?colldiag for more
  ## tips.
  result <- NULL
  if (center) 
    add.intercept <- FALSE
  if (is.matrix(fit) || is.data.frame(fit)) {
    X <- as.matrix(fit)
    nms <- colnames(fit)
  }
  else if (class(fit) == "mer") {
    nms <- names(fixef(fit))
    X <- fit@X
    if (any(grepl("(Intercept)", nms))) {
      add.intercept <- FALSE
    }
  }
  X <- X[!is.na(apply(X, 1, all)), ]
  
  if (add.intercept) {
    X <- cbind(1, X)
    colnames(X)[1] <- "(Intercept)"
  }
  X <- scale(X, scale = scale, center = center)
  
  svdX <- svd(X)
  svdX$d
  condindx <- max(svdX$d)/svdX$d
  dim(condindx) <- c(length(condindx), 1)
  
  Phi = svdX$v %*% diag(1/svdX$d)
  Phi <- t(Phi^2)
  pi <- prop.table(Phi, 2)
  colnames(condindx) <- "cond.index"
  if (!is.null(nms)) {
    rownames(condindx) <- nms
    colnames(pi) <- nms
    rownames(pi) <- nms
  } else {
    rownames(condindx) <- 1:length(condindx)
    colnames(pi) <- 1:ncol(pi)
    rownames(pi) <- 1:nrow(pi)
  }         
  
  result <- data.frame(cbind(condindx, pi))
  zapsmall(result)
}

maxcorr.mer <- function (fit,
                         exclude.intercept = TRUE) {
  so <- summary(fit)
  corF <- so@vcov@factors$correlation
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0 & exclude.intercept) {
    corF <- corF[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  corF[!lower.tri(corF)] <- 0
  maxCor <- max(corF)
  minCor <- min(corF)
  if (abs(maxCor) > abs(minCor)) {
    zapsmall(maxCor)
  } else {
    zapsmall(minCor)
  }
}


vif.mer(lme_3)
kappa.mer(lme_3)

# look at correlation matrix
cor(home_sample_3_long[,c("age","STAI_trait", "pain_cat", "mindfulness", "weight", "cortisol_serum")])


# normality of residuals
describe(residuals(lme_3)) # skew and kurtosis
hist(residuals(lme_3), breaks = 30) # residual histogram
shapiro.test(x=residuals(lme_3)) # not sig
# qq plot
qqmath(lme_3, id=0.05)

#calculate measures of the change in the covariance matrices for the fixed effects based on the deletion of an observation, or group of observations
covratio(lme_3)

#calculate measures of the change in the fixed effects estimates based on the deletetion of an observation, or group of observations
# mdffits(lme_3)

