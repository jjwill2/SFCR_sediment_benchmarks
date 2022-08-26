################################################################################
# Logistic regression between WetSF2.5 and FSBI
# SME50 & SME75 calculation
# Jason Williams, IDEQ Lewiston Office
# last update: 8/5/2022
################################################################################

library(dplyr)
library(ggplot2)

# read in data------------------------------------------------------------------
# Idaho BURP data (all, reference + non-reference)
burp_data <-read.csv("./r_inputs/BURP_data.csv", header = TRUE)
str(burp_data)


# site counts-------------------------------------------------------------------

reg_site_counts <-
  burp_data %>%
  filter(STR_ORDR %in% c("1", "2", "3", "4")) %>%
  filter(!is.na(FSBI)) %>%
  filter(!is.na(wet_fines_2.5_mm)) %>%
  select(BURPID, STR_ORDR, SITECLASS, WBAG.3.Score, wet_fines_2.5_mm, FSBI) %>%
  group_by(SITECLASS, STR_ORDR) %>%
  summarise(site_count = n())

# logistic regression model function--------------------------------------------
# to be applied by siteclass/order combo

logreg_fun <-function(order, siteclass, FSBI_ref) {
  
 
# filter
for_reg <-
  burp_data %>%
  filter(STR_ORDR == order) %>%
  filter(SITECLASS == siteclass) %>%
  filter(!is.na(FSBI)) %>%
  filter(!is.na(wet_fines_2.5_mm)) %>%
  select(BURPID, WBAG.3.Score, wet_fines_2.5_mm, FSBI) %>%
  mutate(FSBI_ref_exceeded = ifelse(FSBI < FSBI_ref, 1, 0))

str(for_reg)

# model
logreg_model <-glm(FSBI_ref_exceeded~wet_fines_2.5_mm, data = for_reg, 
                      family = binomial())
summary(logreg_model)


# chi square test
model_Chi <-logreg_model$null.deviance - logreg_model$deviance
chi_df <-logreg_model$df.null - logreg_model$df.residual
chisq_prob <-1-pchisq(model_Chi, chi_df)

# R2 (Home and Lemeshow)
r2.hl <-model_Chi/logreg_model$null.deviance 
r2.hl

# R2 (Cox and Snell)
r2.cs <-1-exp( -(logreg_model$null.deviance = logreg_model$deviance) / length(logreg_model$fitted.values))
r2.cs

# R2 Naglerke
r2.n <-r2.cs / (1 - ( exp (-(logreg_model$null.deviance / length(logreg_model$fitted.values)))))
r2.n

# odds ratio
exp(logreg_model$coefficients)
exp(confint(logreg_model))

# wet fines at 75% probability of exceeding reference
linedata <-data.frame(wet_fines_2.5_mm = seq(from = 0, to = 100, by = 1))
linedata$predicted = predicted_prob = predict(logreg_model, linedata, type = "response")

percentile50 <-
  linedata %>%
  filter(predicted >=0.50 & predicted <=0.53)

percentile <-
  linedata %>%
  filter(predicted >=0.74 & predicted <=0.76)

# print outputs

print("logistic regression model summary")
print(summary(logreg_model))
print("NOTE: use Mode Chi, not function output summary null & residual deviance")

print("Model Chi")
print(model_Chi)

print("Chi square prob")
print(chisq_prob)

print("Home and Lemeshow R2")
print(r2.hl)

print("wet fines at 51% probability of exceeding reference FSBI")
print(percentile50)

print("wet fines at 75% probability of exceeding reference FSBI")
print(percentile)

# # plot logistic binary data and predicted probability of exceedance

logreg_plot <-
  ggplot() +
  geom_point(data = for_reg, aes(x = wet_fines_2.5_mm, y = FSBI_ref_exceeded)) +
  geom_line(data = linedata, aes(x = wet_fines_2.5_mm, y = predicted_prob)) +
  theme_bw() +
  labs(x = "wetted fines < 2.5 mm", y = "FSBI < reference") +
  ggtitle(paste(siteclass, "", "order ", order))

print(logreg_plot)

} #end logreg fun


# apply function----------------------------------------------------------------

logreg_fun(3, "Foothills", 50)



# plot model & data used to create it
for_mtns_1st_reg_plot <-for_mtns_1st_reg
for_mtns_1st_reg_plot$predicted_prob = predict(mtns_1st_logreg, mtns_1st_logreg_predicted,
                                              type = "response")

linedata <-data.frame(wet_fines_2.5_mm = seq(from = 0, to = 100, by = 1))
linedata$predicted = predicted_prob = predict(mtns_1st_logreg, linedata, type = "response")

plot(predicted~wet_fines_2.5_mm, data = mtns_1st_logreg_predicted)

# plot logistic binary data and predicted probability of exceedance
  
ggplot() +
  geom_point(data = for_mtns_1st_reg, aes(x = wet_fines_2.5_mm, y = FSBI_ref_exceeded)) +
  geom_line(data = linedata, aes(x = wet_fines_2.5_mm, y = predicted_prob)) +
  theme_bw() +
  labs(x = "wetted fines < 2.5 mm", y = "FSBI < reference") 


