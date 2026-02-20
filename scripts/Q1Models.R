# Linear Mixed Models

subsite_table <- readRDS("output/diversitymetrics_table_rare.rds")

library(MASS)
library(lme4)
library(car)
library(multcomp)
library(emmeans)
library(broom)
library(effects)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)

# Tree Density ------------------------------------------------------------

#runnning a linear mixed model to test the effect of green space type on tree abundance, including logged sample area as a covariate because abudance is not rarefied 
modelgst_den0 <- glmer.nb(tree_density_plant_ha ~ GreenSpace + (1|Plot), data= subsite_table)
print(modelgst_den0) #Random effect is zero

#boundary (singular) fit: see help('isSingular') - the random effect variance is near-zero
VarCorr(modelgst_den0) #near-zero, so random effect is not contributing meaningfully
lme4::isSingular(modelgst_den0, tol = 1e-4) 
#Random effect is singular, so remove

#Diagnostics
plot(modelgst_den0, which = 1)
shapiro.test(resid(modelgst_den0))
hist(residuals(modelgst_den0))
qqnorm(residuals(modelgst_den0))
qqline(residuals(modelgst_den0))

#check effect of green space type on abundance using model outputs
Anova(modelgst_den0)
#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.den0 <- broom::tidy(modelgst_den0)
tbl.den0

#tukey test to check pairwise comparisons of abundances of each green space type
Pairwise.den0 <- emmeans(modelgst_den0, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.den0$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.den0) 

#new model - removed random effect 
modelgst_den <- glm.nb(tree_density_plant_ha ~ GreenSpace, data = subsite_table)
print(modelgst_den)
summary(modelgst_den)
confint(modelgst_den, level = 0.95)

# MODEL DIAGNOSTICS
plot(modelgst_den, which = 1)
shapiro.test(resid(modelgst_den))
hist(residuals(modelgst_den))
qqnorm(residuals(modelgst_den))
qqline(residuals(modelgst_den))

#check effect of green space type on tree density using model outputs
Anova(modelgst_den)
#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.den <- broom::tidy(modelgst_den)
tbl.den
plot(allEffects(modelgst_den))
write.csv(tbl.den, "Output/sum_modelgst_abu.csv")

#tukey test to check pairwise comparisons of abundances of each green space type
Pairwise.den <- emmeans(modelgst_den, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.den$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.den) 

# Abundance ---------------------------------------------------------------

#runnning a linear mixed model to test the effect of green space type on tree abundance, including logged sample area as a covariate because abudance is not rarefied 
modelgst_abu0 <- glmer(tree_count ~ GreenSpace + log(area_ha) + (1|Plot), data= subsite_table, family = "negative.binomial"(theta = 1))
print(modelgst_abu0)


#boundary (singular) fit: see help('isSingular') - the random effect variance is near-zero
VarCorr(modelgst_abu0) #near-zero, so random effect is not contributing meaningfully
lme4::isSingular(modelgst_abu0, tol = 1e-4) #Random effect is singular

#Diagnostics
plot(modelgst_abu0, which = 1)
shapiro.test(resid(modelgst_abu0))
hist(residuals(modelgst_abu0))
qqnorm(residuals(modelgst_abu0))
qqline(residuals(modelgst_abu0))

#check effect of green space type on abundance using model outputs
Anova(modelgst_abu0)
#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.abu0 <- broom::tidy(modelgst_abu0)
tbl.abu0

#tukey test to check pairwise comparisons of abundances of each green space type
Pairwise.abu0 <- emmeans(modelgst_abu0, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.abu0$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.abu0) 

#new model - removed random effect 
modelgst_abu <- glm.nb(tree_count ~ GreenSpace + log(area_ha), data = subsite_table)
print(modelgst_abu)
summary(modelgst_abu)
confint(modelgst_abu, level = 0.95)

# Abundance MODEL DIAGNOSTICS
plot(modelgst_abu, which = 1)
shapiro.test(resid(modelgst_abu))
hist(residuals(modelgst_abu))
qqnorm(residuals(modelgst_abu))
qqline(residuals(modelgst_abu))

#check effect of green space type on abundance using model outputs
Anova(modelgst_abu)
#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.abu <- broom::tidy(modelgst_abu)
tbl.abu
plot(allEffects(modelgst_abu))
write.csv(tbl.abu, "Output/sum_modelgst_abu.csv")

#tukey test to check pairwise comparisons of abundances of each green space type
Pairwise.abu <- emmeans(modelgst_abu, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.abu$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.abu) 

#####NOTE
## The pairwise comp. and anova are v diff with random effects and without. 


# Species Richness --------------------------------------------------------

#runnning a linear mixed model to test the effect of green space type on species richness
#do not need to include area as a covariate because we already calculated rarefied diversity metrics based on sample completeness
modelgst_sr0 <- glmer(Cov_rich ~ GreenSpace + (1|Plot), data= subsite_table, family = "negative.binomial"(theta = 1))
print(modelgst_sr0)
#boundary (singular) fit: see help('isSingular') - the random effect variance is near-zero
VarCorr(modelgst_sr0) #near-zero, so random effect is not contributing meaningfully
lme4::isSingular(modelgst_sr0, tol = 1e-4) #Random effect is singular

plot(modelgst_sr0, which = 1)
shapiro.test(resid(modelgst_sr0))
hist(residuals(modelgst_sr0))
qqnorm(residuals(modelgst_sr0))
qqline(residuals(modelgst_sr0))

#new model - removed random effect 
modelgst_sr <- glm.nb(Cov_rich ~ GreenSpace, data = subsite_table)
print(modelgst_sr)
summary(modelgst_sr)
confint(modelgst_sr, level = 0.95)

# SR MODEL DIAGNOSTICS
plot(modelgst_sr, which = 1)
shapiro.test(resid(modelgst_sr))
hist(residuals(modelgst_sr))
qqnorm(residuals(modelgst_sr))
qqline(residuals(modelgst_sr))

#check effect of green space type on species richness using model outputs
Anova(modelgst_sr)
#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.sr <- broom::tidy(modelgst_sr)
tbl.sr
plot(allEffects(modelgst_sr))
write.csv(tbl.sr, "Output/sum_modelgst_sr.csv")

#tukey test to check pairwise comparisons of species richness of each green space type
Pairwise.sr <- emmeans(modelgst_sr, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.sr$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.sr) 


# Effective Number of Species ---------------------------------------------

#runnning a linear mixed model to test the effect of green space type on ENS
#do not need to include area as a covariate because we already calculated rarefied diversity metrics based on sample completeness
modelgst_ens <- lmer(Cov_q1 ~ GreenSpace + (1|Plot), data= subsite_table)
summary(modelgst_ens) #Random effect explains some variance, will include

# ENS MODEL DIAGNOSTICS
qqnorm(residuals(modelgst_ens))
qqline(residuals(modelgst_ens))
hist(residuals(modelgst_ens))
plot(residuals(modelgst_ens))
shapiro.test(resid(modelgst_ens))

modelgst_ens_glm <- glm.nb(Cov_q1 ~ GreenSpace, data = subsite_table)

qqnorm(residuals(modelgst_ens_glm))
qqline(residuals(modelgst_ens_glm))
hist(residuals(modelgst_ens_glm))
plot(residuals(modelgst_ens_glm))
shapiro.test(resid(modelgst_ens_glm))


#Residuals are not normal, so try Gamma GLMM with log link
modelgst_ens_gamma <- glmer(Cov_q1 ~ GreenSpace + (1 | Plot), data = subsite_table, family = Gamma(link = "log"))
summary(modelgst_ens_gamma) #Random effect is small

#Residuals
resid_gamma <- residuals(modelgst_ens_gamma, type = "pearson")
qqnorm(resid_gamma); qqline(resid_gamma)
hist(resid_gamma)
plot(resid_gamma)

#Compare with gamma glm - no random effect
modelgst_ens_gamma_glm <- glm(Cov_q1 ~ GreenSpace, data = subsite_table, family = Gamma(link = "log"))
summary(modelgst_ens_gamma_glm)

qqnorm(residuals(modelgst_ens_gamma_glm))
qqline(residuals(modelgst_ens_gamma_glm))
hist(residuals(modelgst_ens_gamma_glm))
plot(residuals(modelgst_ens_gamma_glm))
shapiro.test(resid(modelgst_ens_gamma_glm))


AIC(modelgst_ens_gamma, modelgst_ens_gamma_glm) #gamma glmer has margianlly lower AIC, so go with this model

#check effect of green space type on species richness using model outputs
Anova(modelgst_ens_gamma)
#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.ens <- tidy(modelgst_ens_gamma)
tbl.ens
plot(allEffects(modelgst_ens_gamma))
write.csv(tbl.ens, "Output/sum_modelgst_ens.csv")

#tukey test to check pairwise comparisons of species richness of each green space type
Pairwise.ens <- emmeans(modelgst_ens_gamma, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.ens$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.ens) 


# Evenness ----------------------------------------------------------------

#runnning a linear mixed model to test the effect of green space type on evenness
modelgst_eve <- lmer(Cov_pie ~ GreenSpace + (1|Plot), data= subsite_table)
summary(modelgst_eve)

# Evenness MODEL DIAGNOSTICS
qqnorm(residuals(modelgst_eve))
qqline(residuals(modelgst_eve))
hist(residuals(modelgst_eve))
plot(residuals(modelgst_eve))
shapiro.test(resid(modelgst_eve)) #Residuals are not normal

#Try beta regression, mixed effects
# Small constant to shrink 0s and 1s inward, because beta regression only handles values between 0 and 1. 
eps <- 0.001
subsite_table$Cov_pie_adj <- subsite_table$Cov_pie
# Shrink 0s up to eps
subsite_table$Cov_pie_adj[subsite_table$Cov_pie_adj == 0] <- eps

modelgst_eve_beta <- glmmTMB(Cov_pie_adj ~ GreenSpace + (1 | Plot), family = beta_family(link = "logit"), data = subsite_table)
summary(modelgst_eve_beta) #Random effects explain variance, will keep

#Check model fit
# Simulate residuals from the beta GLMM
sim_res <- simulateResiduals(fittedModel = modelgst_eve_beta)
# Plot residual diagnostics
plot(sim_res)
# Perform formal tests for uniformity, dispersion, zero-inflation etc.
testUniformity(sim_res)
testDispersion(sim_res)
#NOT A GOOD FIT

#Also tried Logit Transformation + Linear Mixed Model (LMM), was not a good fit

# Fit zero inflated beta regression
modelgst_eve_beta_01 <- glmmTMB(Cov_pie ~ GreenSpace + (1 | Plot), 
                             family = beta_family(link = "logit"),
                             ziformula = ~1,  
                             data = subsite_table
)

summary(modelgst_eve_beta_01)

# Diagnostics with DHARMa
sim_res <- simulateResiduals(modelgst_eve_beta_01)
plot(sim_res)
testUniformity(sim_res)
testDispersion(sim_res)

#check effect of green space type on evenness using model outputs
Anova(modelgst_eve_beta_01)
#Green space statistically significant portion of the variation
summary(modelgst_eve_beta_01)$coefficients$cond
VarCorr(modelgst_eve_beta_01)$cond

#tukey test to check pairwise comparisons of evenness of each green space type
Pairwise.eve <- emmeans(modelgst_eve_beta_01, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.eve$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.eve) 


