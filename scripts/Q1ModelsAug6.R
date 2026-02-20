#NEW NEW Q1 MODELS AUGUST 6th 2025

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

modelgst_den <- glmer.nb(tree_count ~ GreenSpace + offset(log(plantarea_ha)) + (1|Plot), data= subsite_table)
print(modelgst_den) 
summary(modelgst_den)

#Diagnostics
plot(modelgst_den, which = 1)
shapiro.test(resid(modelgst_den))
hist(residuals(modelgst_den))
qqnorm(residuals(modelgst_den))
qqline(residuals(modelgst_den))

#Just ran this to compare with below
#check effect of green space type on abundance using model outputs
Anova(modelgst_den)
#Green space statistically significant portion of the variation
tbl.den <- broom::tidy(modelgst_den)
tbl.den

#tukey test to check pairwise comparisons of abundances of each green space type
Pairwise.den <- emmeans(modelgst_den, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.den$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.den) 

# Species Richness --------------------------------------------------------

#Gamma with log link
modelgst_sr <- glmer(Cov_rich ~ GreenSpace + (1 | Plot), data = subsite_table, family = Gamma(link = "log"))
summary(modelgst_sr) 

#Residuals
resid_gamma <- residuals(modelgst_sr, type = "pearson")
qqnorm(resid_gamma); qqline(resid_gamma)
hist(resid_gamma)
plot(resid_gamma)

#check effect of green space type on species richness using model outputs
Anova(modelgst_sr)
#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.sr <- tidy(modelgst_sr)
tbl.sr

#tukey test to check pairwise comparisons of species richness of each green space type
Pairwise.sr <- emmeans(modelgst_sr, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.sr$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.sr) 

# Effective Number of Species ---------------------------------------------

#Gamma with log link
modelgst_ens <- glmer(Cov_q1 ~ GreenSpace + (1 | Plot), data = subsite_table, family = Gamma(link = "log"))
summary(modelgst_ens) 

#Residuals
resid_gamma <- residuals(modelgst_ens, type = "pearson")
qqnorm(resid_gamma); qqline(resid_gamma)
hist(resid_gamma)
plot(resid_gamma)

#check effect of green space type on ENS using model outputs
Anova(modelgst_ens)
#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.ens <- tidy(modelgst_ens)
tbl.ens

#tukey test to check pairwise comparisons of ENS of each green space type
Pairwise.ens <- emmeans(modelgst_ens, pairwise ~ GreenSpace, type = "response")
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

#Green space statistically significant portion of the variation
#push out model summary as csv file for plotting
tbl.eve <- tidy(modelgst_eve_beta_01)
tbl.eve

#tukey test to check pairwise comparisons of evenness of each green space type
Pairwise.eve <- emmeans(modelgst_eve_beta_01, pairwise ~ GreenSpace, type = "response")
# Summary of pairwise comparisons
summary(Pairwise.eve$contrasts)
#summarizing pairwise differences across green space types
summary(Pairwise.eve) 


###Modelling tree density:

newdat <- expand.grid(GreenSpace = unique(subsite_table$GreenSpace),
                      plantarea_ha = 1)

# Get predictions
preds <- predict(modelgst_den,
                 newdata = newdat,
                 re.form = NA,   # fixed effects only
                 type = "response",
                 se.fit = TRUE)

pred_df <- newdat %>%
  mutate(mean = preds$fit,
         se = preds$se.fit)

# Fixed order for GreenSpace types
ordered_levels <- c("Commercial", "Vacant Lot", "Park",
                    "Institutional", "Public Right-Of-Way", "Residential")

pred_df <- pred_df %>%
  mutate(GreenSpace = factor(GreenSpace, levels = ordered_levels))

# Plot
ggplot(pred_df, aes(x = GreenSpace, y = mean, color = GreenSpace)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.5, linewidth = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(x = "Green Space Type",
       y = "Tree Density (per ha)",
       title = "Model-predicted Tree Density by Green Space Type")
