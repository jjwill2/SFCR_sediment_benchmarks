################################################################################
# Part 1: exploratory data analysis/plots
# Part 2: WetSF2.5 and FSBI percentile-based reference benchmarks
# Jason Williams, IDEQ Lewiston Office
# last update: 8/5/2022
################################################################################


library(dplyr)
library(ggplot2)
library(reshape2)

# read in data------------------------------------------------------------------

# Idaho BURP data (all, reference + non-reference)
burp_data <-read.csv("./r_inputs/BURP_data.csv", header = TRUE)
str(burp_data)

################################################################################
# Part 1: Exploratory Analysis
################################################################################

# count of BURP sites by WBAG 3 siteclass and reference status------------------

burp_data_formatted <-
  burp_data %>%
  select(BURPID, REFERENCE, SITECLASS, STR_ORDR, FSBI, wet_fines_2.5_mm, SMI2) %>%
  melt(id.vars = c("BURPID", "SITECLASS", "REFERENCE", "STR_ORDR"), 
       variable.name = "parameter", value.name = "value") %>%
  filter(!is.na(value))

str(burp_data_formatted)

burp_site_counts <-
  burp_data_formatted %>%
  group_by(parameter, SITECLASS, REFERENCE) %>%
  summarise(site_count = n())

# also by order
burp_site_counts_order <-
  burp_data_formatted %>%
  group_by(parameter, SITECLASS, REFERENCE, STR_ORDR) %>%
  summarise(site_count = n())

# all reference data across siteclasses by order
burp_data_formatted %>%
  filter(REFERENCE == "Reference") %>%
  filter(STR_ORDR %in% c("1", "2", "3", "4")) %>%
  filter(!is.na(value)) %>%
  group_by(parameter, STR_ORDR) %>%
  summarise(site_count = n())

# all reference data across order by siteclass
burp_data_formatted %>%
  filter(REFERENCE == "Reference") %>%
  filter(STR_ORDR %in% c("1", "2", "3", "4")) %>%
  filter(!is.na(value)) %>%
  group_by(parameter, SITECLASS) %>%
  summarise(site_count = n())


# correlation matrix------------------------------------------------------------

burp_for_corr <-
  burp_data %>%
  filter(STR_ORDR != "0L") %>%
  mutate(order_numeric = as.numeric(STR_ORDR), 
         elevation_numeric = as.numeric(ELEVATION),
         catchment_numeric = as.numeric(CATCHMENT),
         FSBI_numeric = as.numeric(FSBI),
         SMI2_numeric = as.numeric(SMI2)) %>%
  select(REFERENCE, SITECLASS, order_numeric, GRADIENT, elevation_numeric, 
         catchment_numeric, 
         FSBI_numeric, wet_fines_2.5_mm, all_fines_2.5_mm, SMI, WBAG.2.Score,
         SMI2_numeric, WBAG.3.Score, riffle_run_pct, riffle_pct) # NA warnings ok

str(burp_for_corr)

# correlation matrix all data
for_corr1 <-burp_for_corr %>% select(-REFERENCE, -SITECLASS)
corr_matrix <-cor(for_corr1, use = "pairwise.complete.obs", 
                  method = "pearson")
round(corr_matrix, 2)

# correlation matrix, reference data only
for_corr2 <-
  burp_for_corr %>% 
  filter(REFERENCE == "Reference") %>%
  select(-REFERENCE, -SITECLASS)
ref_corr_matrix <-cor(for_corr2, use = "pairwise.complete.obs", method = "pearson")
round(ref_corr_matrix, 2)


# correlation matrix, mountains reference data only
for_corr3 <-
  burp_for_corr %>% 
  filter(REFERENCE == "Reference") %>%
  filter(SITECLASS == "Mountains") %>%
  select(-REFERENCE, -SITECLASS)
mtn_corr_matrix <-cor(for_corr3, use = "pairwise.complete.obs", method = "pearson")
round(mtn_corr_matrix, 2)

# correlation matrix, PPBV reference data only
for_corr4 <-
  burp_for_corr %>% 
  filter(REFERENCE == "Reference") %>%
  filter(SITECLASS == "PPBV") %>%
  select(-REFERENCE, -SITECLASS)
ppbv_corr_matrix <-cor(for_corr4, use = "pairwise.complete.obs", method = "pearson")
round(ppbv_corr_matrix, 2)

# correlation matrix, PPBV reference data only
for_corr5 <-
  burp_for_corr %>% 
  filter(REFERENCE == "Reference") %>%
  filter(SITECLASS == "Foothills") %>%
  select(-REFERENCE, -SITECLASS)
foot_corr_matrix <-cor(for_corr5, use = "pairwise.complete.obs", method = "pearson")
round(ppbv_corr_matrix, 2)


# fines exploratory plots-------------------------------------------------------


# BURP wet fines < 2.5 mm box plot
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE %in% c("Reference", "Stress")) %>%
  ggplot(aes(x = REFERENCE, y = wet_fines_2.5_mm)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey") +
  facet_grid(~SITECLASS)

# wet fines by order, reference only, faceted by siteclass
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = STR_ORDR, y = wet_fines_2.5_mm)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey") +
  facet_grid(~SITECLASS) +
  theme_bw() +
  labs(x = "stream order", y = "wetted width % fines < 2.5 mm")

# wet fines by order, reference only, NOT faceted by siteclass
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = STR_ORDR, y = wet_fines_2.5_mm)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey") +
  theme_bw() +
  #facet_grid(~SITECLASS) +
  labs(x = "stream order", y = "wetted width % fines < 2.5 mm")

# wet fines response to gradient
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = GRADIENT, y = wet_fines_2.5_mm)) +
  geom_point() +
  facet_grid(SITECLASS~REFERENCE)

# wet fines response to gradient by order, reference only
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE =="Reference") %>%
  ggplot(aes(x = GRADIENT, y = wet_fines_2.5_mm)) +
  geom_point() +
  facet_grid(SITECLASS~STR_ORDR)

# wet fines vs % riffle/run
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE =="Reference") %>%
  ggplot(aes(x = riffle_run_pct, y = wet_fines_2.5_mm)) +
  geom_point()

# wet fines vs % riffle
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE =="Reference") %>%
  ggplot(aes(x = riffle_pct, y = wet_fines_2.5_mm)) +
  geom_point()

# BURP wet fines < 2.5 mm box plot by habitat type, reference only
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  mutate(hab_type = ifelse(riffle_run_pct < 35, "pool-glide", 
                           ifelse(riffle_run_pct >= 65, "riffle/run", "other"))) %>%
  ggplot(aes(x = hab_type, y = wet_fines_2.5_mm)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey") +
  facet_grid(~SITECLASS)

# FSBI exploratory plots--------------------------------------------------------

# FSBI box plot
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE %in% c("Reference", "Stress")) %>%
  ggplot(aes(x = REFERENCE, y = FSBI)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey") +
  facet_grid(~SITECLASS)

# FSBI box plot, reference only
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = SITECLASS, y = FSBI)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey") +
  theme_bw()

# FSBI by order, reference only
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = STR_ORDR, y = FSBI)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey") +
  facet_grid(~SITECLASS) 

# FSBI response to wet fines (reference only)
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = wet_fines_2.5_mm, y = FSBI)) +
  geom_point() +
  geom_smooth() +
  facet_grid(~SITECLASS)

# FSBI response to wet fines by order (reference only)
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = wet_fines_2.5_mm, y = FSBI)) +
  geom_point() + 
  facet_grid(SITECLASS~STR_ORDR)

# SMI2 vs FSBI, reference only
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = SMI2, y = FSBI)) +
  geom_point() +
  facet_grid(~SITECLASS) +
  theme_bw() 

# % riffle vs FSBI, reference only, by siteclass
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = riffle_pct, y = FSBI)) +
  geom_point() +
  facet_grid(~SITECLASS)

# % riffle/run vs FSBI, reference only, by siteclass
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  ggplot(aes(x = riffle_run_pct, y = FSBI)) +
  geom_point() +
  facet_grid(~SITECLASS)

# % pool vs FSBI, reference only, by siteclass
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  mutate(pool_pct = (pool / total) * 100) %>%
  ggplot(aes(x = pool_pct, y = FSBI)) +
  geom_point() +
  facet_grid(~SITECLASS)

# FSBI riffle/run vs glide/pool, reference only by siteclass
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  filter(!is.na(riffle_run_pct)) %>%
  mutate(hab_type = ifelse(riffle_run_pct < 35, "pool-glide", 
                           ifelse(riffle_run_pct >= 35, "riffle/run", "other"))) %>%
  ggplot(aes(x = hab_type, y = FSBI)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey") +
  facet_grid(~SITECLASS) +
  theme_bw() +
  labs(x = "reach habitat type", y = "FSBI")

# FSBI riffle/run vs glide/pool, reference only all siteclasses
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(!is.na(riffle_run_pct)) %>%
  filter(REFERENCE == "Reference") %>%
  mutate(hab_type = ifelse(riffle_run_pct < 35, "pool-glide", 
                           ifelse(riffle_run_pct >= 35, "riffle/run", "other"))) %>%
  ggplot(aes(x = hab_type, y = FSBI)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey")

# FSBI riffle/run vs not riffle/run
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  filter(!is.na(riffle_run_pct)) %>%
  mutate(hab_type = ifelse(riffle_run_pct >= 35, "riffle/run", "other")) %>%
  ggplot(aes(x = hab_type, y = FSBI)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey")

# FSBI riffle/run vs not riffle/run
burp_data %>%
  filter(STR_ORDR < 5) %>%
  filter(REFERENCE == "Reference") %>%
  filter(!is.na(riffle_run_pct)) %>%
  mutate(pool_pct = (pool / total) * 100) %>%
  mutate(hab_type = ifelse(riffle_pct > 50, "mostly riffle", 
                           ifelse(pool_pct > 50, "mostly pool", "mix"))) %>%
  ggplot(aes(x = hab_type, y = FSBI)) +
  geom_boxplot() +
  geom_point(shape = 1, color = "darkgrey")


################################################################################
# Part 2 Assessment Threshold Calculations
################################################################################

# Reference percentiles, % fines------------------------------------------------

# all reference data, wet % fines < 2.5 mm
ref <-
  burp_data %>%
  filter(REFERENCE == "Reference")

# 75th percentile, all data across orders
quantile(ref$wet_fines_2.5_mm, probs=0.75, na.rm = TRUE) 

# 90th percentile
quantile(ref$wet_fines_2.5_mm, probs=0.9, na.rm = TRUE) 

# 75th percentile, 1st order 
ref_1st <-ref %>% filter(STR_ORDR == "1")
quantile(ref_1st$wet_fines_2.5_mm, probs=0.75, na.rm = TRUE) 

# 75th percentile, 2nd order 
ref_2nd <-ref %>% filter(STR_ORDR == "2")
quantile(ref_2nd$wet_fines_2.5_mm, probs=0.75, na.rm = TRUE) 

# 75th percentile, 3rd order
ref_3rd <-ref %>% filter(STR_ORDR == "3")
quantile(ref_3rd$wet_fines_2.5_mm, probs=0.75, na.rm = TRUE) 

# 75th percentile, 4th order
ref_4th <-ref %>% filter(STR_ORDR == "4")
quantile(ref_4th$wet_fines_2.5_mm, probs=0.75, na.rm = TRUE) 

# Reference percentiles, FSBI---------------------------------------------------
# 25th percentile because FSBI decreases with increasing % fines

# mtn siteclass 25th percentile of reference
mtn_ref <-
  burp_data %>%
  filter(REFERENCE == "Reference") %>%
  filter(SITECLASS == "Mountains")
quantile(mtn_ref$FSBI, probs=0.25, na.rm = TRUE) 

# foothills siteclass 75th percentile of reference
foot_ref <-
  burp_data %>%
  filter(REFERENCE == "Reference") %>%
  filter(SITECLASS == "Foothills")
quantile(foot_ref$FSBI, probs=0.25, na.rm = TRUE) 

# PPBV siteclass 75th percentile of reference
ppbv_ref <-
  burp_data %>%
  filter(REFERENCE == "Reference") %>%
  filter(SITECLASS == "PPBV")
quantile(ppbv_ref$FSBI, probs=0.25, na.rm = TRUE) 
