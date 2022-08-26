################################################################################
# Evaluates sediment assessment protocol performance as described in 
# APpendix A section A3
# Jason Williams, IDEQ Lewiston Office
# last update: 8/5/2022
################################################################################


library(dplyr)
library(ggplot2)

# read in data------------------------------------------------------------------
# Idaho BURP data (all, reference + non-reference)
burp_data <-read.csv("./r_inputs/BURP_data.csv", header = TRUE)
str(burp_data)

cat4a <-read.csv("./r_inputs/cat4a_2022IR.csv", header = TRUE)
str(cat4a)

bench <-read.csv("./r_inputs/benchmarks_for_assessment.csv", header = TRUE)
str(bench)


# format TMDL data--------------------------------------------------------------

levels(unique(as.factor(cat4a$CAUSE_NAME)))

cat4a_sed <-
  cat4a %>%
  filter(CAUSE_NAME %in% c("SEDIMENTATION/SILTATION", 
                           "TOTAL SUSPENDED SOLIDS (TSS)"))
  

# evaluate if benchmarks exceeded-----------------------------------------------


assessment_data_formatted <-
  burp_data %>%
  filter(STR_ORDR %in% c("1", "2", "3", "4")) %>%
  select(BURPID, YEAR, SITECLASS, WBSEGMENT, STR_ORDR, FSBI, wet_fines_2.5_mm, SMI2, 
         WBAG.3.Score) %>%
  mutate(order = as.numeric(STR_ORDR)) %>%
  merge(bench, by.x = c("SITECLASS", "order"), by.y = c("Siteclass", "Order"), 
        all.x = TRUE) %>%
  
  # has sed TMDL?
  mutate(sed_TMDL = ifelse(WBSEGMENT %in% cat4a_sed$ID305B, "Y", "N")) %>%
  
  # benchmarks exceeded?
  mutate(wetSF2.5_ref_exceeded = ifelse(wet_fines_2.5_mm > wetSF2.5_ref_benchmark,
         1, 0)) %>%
  mutate(FSBI_ref_exceeded = ifelse(FSBI > FSBI_ref_benchmark, 1, 0)) %>%
  mutate(SME75_exceeded = ifelse(wet_fines_2.5_mm > SME75_benchmark, 1, 0)) %>%
  mutate(SMI2_fail = ifelse(SMI2 > SMI2_WBAG3_score2, 0, 1)) %>%
  
  # ref exceeded?
  mutate(ref_score = wetSF2.5_ref_exceeded + FSBI_ref_exceeded) %>%
  
  # total sed score 
  mutate(sed_score = ref_score + SME75_exceeded + SMI2_fail)

str(assessment_data_formatted)


# False positive test1 ---------------------------------------------------------
# false positive if BURP & SMI2 pass, but sed impairment predicted

false_pos_test <-
  assessment_data_formatted %>%
  filter(!is.na(WBAG.3.Score)) %>%
  filter(!is.na(FSBI)) %>%
  filter(!is.na(wet_fines_2.5_mm)) %>%
  filter(!is.na(SMI2)) %>%
  mutate(sed_impair_predicted = ifelse(sed_score >=3, "sediment impairment predicted",
                                       "sediment impairment not predicted")) %>%
  mutate(status = ifelse(WBAG.3.Score >= 2 & SMI2_fail == 0, "BURP & SMI2 Pass", 
                         ifelse(WBAG.3.Score <2 & SMI2_fail == 1, "BURP & SMI2 Fail", "other")))

# test 1
false_pos_test %>%
  group_by(status, sed_impair_predicted) %>%
  summarise(count = n())

# test 2
false_pos_test %>%
  group_by(SMI2_fail, sed_impair_predicted) %>%
  summarise(count = n())


# false negative test-----------------------------------------------------------

false_neg_test <-
  assessment_data_formatted %>%
  filter(YEAR >= 2015) %>%
  filter(!is.na(FSBI)) %>%
  filter(!is.na(wet_fines_2.5_mm)) %>%
  filter(!is.na(SMI2)) %>%
  mutate(sed_impair_predicted = ifelse(sed_score >=3, "sediment impairment predicted",
                                       "sediment impairment not predicted")) 

false_neg_test %>%
  group_by(sed_TMDL, sed_impair_predicted) %>%
  summarise(count = n())
  
  