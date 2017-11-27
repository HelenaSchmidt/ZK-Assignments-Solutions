#### Solution ZK Assignment 1 ####
library("psych")
library("car")
library("lm.beta")
library("lsr")

# read in data
home_sample_1 = read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class/master/home_sample_1.csv") 

# copy data
home_sample_11 = home_sample_1
View(home_sample_11)

# data diagnostics and identification of coding errors
describe(home_sample_11)
summary(home_sample_11)
home_sample_11=home_sample_11[!home_sample_11$sex=="3",]
boxplot(home_sample_11$mindfulness)
stem(home_sample_11$mindfulness)
home_sample_11=home_sample_11[!home_sample_11$mindfulness<1,]
home_sample_11[!complete.cases(home_sample_11),] # no missing values found

# histograms
hist(home_sample_11$age, breaks = 30)
hist(home_sample_11$STAI_trait, breaks = 30)
hist(home_sample_11$pain_cat, breaks = 30)
hist(home_sample_11$mindfulness, breaks = 30)
hist(home_sample_11$cortisol_serum, breaks = 30)
hist(home_sample_11$cortisol_saliva, breaks = 30)

# boxplot
boxplot(home_sample_11$age)
boxplot(home_sample_11$STAI_trait)
boxplot(home_sample_11$pain_cat)
boxplot(home_sample_11$cortisol_saliva)

# scatterplots
plot(pain ~ sex, data = home_sample_11)
plot(pain ~ age, data = home_sample_11)
abline(lm(pain~ age, data=home_sample_11))
plot(pain ~ STAI_trait, data = home_sample_11)
abline(lm(pain~ STAI_trait, data=home_sample_11))
plot(pain ~ pain_cat, data = home_sample_11)
abline(lm(pain~ pain_cat, data=home_sample_11))
plot(pain ~ mindfulness, data = home_sample_11)
abline(lm(pain~ mindfulness, data=home_sample_11))
plot(pain ~ cortisol_serum, data = home_sample_11)
abline(lm(pain~ cortisol_serum, data=home_sample_11))
plot(pain ~ cortisol_saliva, data = home_sample_11)
abline(lm(pain~ cortisol_saliva, data=home_sample_11))

# theoretically established model 1
mod1 <- lm(pain ~ sex + age, data = home_sample_11)
summary(mod1)
AIC(mod1)
BIC(mod1)
confint(mod1)
lm.beta(mod1)
standardCoefs(mod1)

# model 2 with psychological and hormonal predictors
mod2 = lm(pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_serum + cortisol_saliva, data = home_sample_11)
summary(mod2)
AIC(mod2)
confint(mod2)
lm.beta(mod2)
standardCoefs(mod2)

# assumption-testing model 1
# multicollinearity
vif(mod1)

# normality of residuals
plot(x=mod1, which = 2) # QQ
describe(residuals(mod1)) # skew and kurtosis
hist(residuals(mod1), breaks = 20) # residual histogram
shapiro.test(x=residuals(mod1)) # sig!

# linearity
pred <- predict(object = mod1) # predicted values against actual values
plot( x = pred, y = home_sample_11$pain,
      xlab = "Fitted Values", ylab = "Observed Values")
plot(x=mod1, which = 1) # predicted values against residuals
residualPlots(mod1)

# homoscedasticity
plot(mod1, which = 3)
ncvTest(mod1) # not sig

# identifiy influential outliers
# Cooks distance
plot(cooks.distance(mod1))
plot(mod1, pch=18)
plot(mod1, which=4)
rev(sort(cooks.distance(mod1)))

# assumption-testing model 2
# multicollinearity
vif(mod2) # OBS! cortisol_serum and cortisol_saliva over 5.
pairs.panels(home_sample_11[,c("pain", "sex", "age", "STAI_trait", "pain_cat", "mindfulness", "cortisol_serum", "cortisol_saliva")], col = "red", lm = T)
# look at correlation between cortisol measurements
cor(home_sample_11$cortisol_serum, home_sample_11$cortisol_saliva) # highly correlated! .87

# compare the two cortisol measurements
mod2_new1 = lm(pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = home_sample_11)
mod2_new2 = lm(pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_saliva, data = home_sample_11)
summary(mod2_new1)
summary(mod2_new2)
AIC(mod2_new1)
AIC(mod2_new2)
BIC(mod2_new1)
BIC(mod2_new2) # decision: choose model new1 with cortisol_saliva if looking at anova and AIC, BUT:
# OBS! But instead of saliva we choose the serum measurement, because it is more used in the theory
# OBS! We avoid overfitting of our model
confint(mod2_new1)
lm.beta(mod2_new1)
standardCoefs(mod2_new1) # include zero for the variables without significant effect in the model!


# test multicollinearity again
vif(mod2_new1)
pairs.panels(home_sample_11[,c("pain", "sex", "age", "STAI_trait", "pain_cat", "mindfulness", "cortisol_serum")], col = "red", lm = T)

# normality of residuals
plot(x=mod2_new1, which = 2) # QQ
describe(residuals(mod2_new1)) # skew and kurtosis
hist(residuals(mod2_new1), breaks = 20) # residual histogram
shapiro.test(x=residuals(mod2_new1)) # not sig

# linearity
pred <- predict(object = mod2_new1) # predicted values against actual values
plot( x = pred, y = home_sample_11$pain,
      xlab = "Fitted Values", ylab = "Observed Values")
plot(x=mod2_new1, which = 1) # predicted values against residuals
residualPlots(mod2_new1)

# homoscedasticity
plot(mod2_new1, which = 3)
ncvTest(mod2_new1) # not sig

# identifiy influential outliers
# Cooks distance
plot(cooks.distance(mod2_new1))
plot(mod2_new1, pch=18)
plot(mod2_new1, which=4)
rev(sort(cooks.distance(mod2_new1)))

# model comparison (nested models)
summary(mod1)
confint(mod1)
lm.beta(mod1)
standardCoefs(mod1)
summary(mod2_new1)
confint(mod2_new1)
lm.beta(mod2_new1)
standardCoefs(mod2_new1)

anova(mod1, mod2_new1)
AIC(mod1)
AIC(mod2_new1)
BIC(mod1)
BIC(mod2_new1)
