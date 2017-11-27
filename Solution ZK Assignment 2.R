#### Solution ZK Assignment 1 ####
library("psych")
library("car")
library("lm.beta")

# clear environment
remove(list = ls())

# read in data
home_sample_1 = read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class/master/home_sample_1.csv") 

# copy data
home_sample_11 = home_sample_1


# exclude cases like in Assignment 1
home_sample_11=home_sample_11[!home_sample_11$sex=="3",]
home_sample_11=home_sample_11[!home_sample_11$mindfulness<1,]

# data diagnostics of new variable weight
boxplot(home_sample_11$weight)
hist(home_sample_11$weight, breaks = 30)
plot(pain ~ weight, data = home_sample_11)

# initial regression for backward regression
mod_initial = lm(formula= pain ~ age + sex + weight + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = home_sample_11)

# assumption-testing
# normality of residuals
plot(x=mod_initial, which = 2) # QQ
describe(residuals(mod_initial)) # skew and kurtosis
hist(residuals(mod_initial), breaks = 20) # residual histogram
shapiro.test(x=residuals(mod_initial)) # not sig.

# linearity
pred <- predict(object = mod_initial ) # predicted values against actual values
plot( x = pred, y = home_sample_11$pain,
      xlab = "Fitted Values", ylab = "Observed Values")
plot(x=mod_initial, which = 1) # predicted values against residuals
residualPlots(mod_initial)

# homoscedasticity
plot(mod_initial, which = 3)
ncvTest(mod_initial) # not sig.

# multicollinearity
vif(mod_initial)
pairs.panels(home_sample_11[,c("pain", "age", "sex", "weight", "STAI_trait", "pain_cat", "mindfulness", "cortisol_serum")], col = "red", lm = T)

# Cooks distance
plot(cooks.distance(mod_initial))
plot(mod_initial, pch=18)
plot(mod_initial, which=4)
rev(sort(cooks.distance(mod_initial)))  # highest rows: 37, 123, 160

# backward regression
mod_back1 = step(mod_initial, direction = "backward")

# backward model
mod_backward <- lm(formula=pain ~ age + pain_cat + mindfulness + cortisol_serum, data=home_sample_11)
summary(mod_backward)
confint(mod_backward)
lm.beta(mod_backward)
standardCoefs(mod_backward)
AIC(mod_backward)
AIC(mod_initial)

# assumption-testing of backward model
# normality of residuals
plot(x=mod_backward, which = 2) # QQ
describe(residuals(mod_backward)) # skew and kurtosis
hist(residuals(mod_backward), breaks = 20) # residual histogram
shapiro.test(x=residuals(mod_backward)) # not sig.

# linearity
pred <- predict(object = mod_backward ) # predicted values against actual values
plot( x = pred, y = home_sample_11$pain,
      xlab = "Fitted Values", ylab = "Observed Values")
plot(x=mod_backward, which = 1) # predicted values against residuals
residualPlots(mod_backward)

# homoscedasticity
plot(mod_initial, which = 3)
ncvTest(mod_backward) # not sig.

# multicollinearity
vif(mod_initial)
pairs.panels(home_sample_11[,c("pain", "age", "sex", "weight", "STAI_trait", "pain_cat", "mindfulness", "cortisol_serum")], col = "red", lm = T)

#4 plots in one mod_backward
par(mfrow=c(2,2))
plot(mod_backward)
par(mfrow=c(1,1))

# Cooks distance
plot(cooks.distance(mod_backward))
plot(mod_backward, pch=18)
plot(mod_backward, which=4)
rev(sort(cooks.distance(mod_backward)))

# theory based model
mod_theorybased <- lm(formula=pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = home_sample_11)

# model comparsion full model and backward model
anova(mod_backward, mod_initial) # no significant difference, but simpler model prefered
summary(mod_inital)
summary(mod_backward)
AIC(mod_initial)
AIC(mod_backward)
BIC(mod_initial)
BIC(mod_backward)

# model comparison
anova(mod_backward, mod_theorybased)
summary(mod_backward)
summary(mod_theorybased)
standardCoefs(mod_backward)
standardCoefs(mod_theorybased)
AIC(mod_backward)
AIC(mod_theorybased)
BIC(mod_backward)
BIC(mod_theorybased)

# Test models on new data
# read in data
home_sample_2 = read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class/master/home_sample_2.csv")

# data diagnostics and identification of coding errors
describe(home_sample_2)
summary(home_sample_2)
boxplot(home_sample_2$mindfulness)
stem(home_sample_2$mindfulness)

# copy data
home_sample21 = home_sample_2

#remove cases with mindufulness values smaller than 1
home_sample_21=home_sample_2[!home_sample_2$mindfulness<1,]
home_sample_21[!complete.cases(home_sample_21),] # no missing values found

summary(home_sample_21)

# predict pain on new data without re-fitting models
pred_backward <- predict(mod_backward, newdata=home_sample_21)
pred_theorybased <- predict(mod_theorybased, newdata=home_sample_21)
summary(pred_backward)
summary(pred_theorybased)

# calculate sum of squared residuals - where is more error?
RSS_backward = sum((home_sample_21$pain - pred_backward)^2)
RSS_theorybased = sum((home_sample_21$pain - pred_theorybased)^2)
RSS_backward
RSS_theorybased
