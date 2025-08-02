## **
# Following Qi's email on 05/03/2023 (subject: [Encrypt] CCSS help). Running Piecewise-exponential regression.
rm(list=ls())
#########################
## Load Phenotype data ##
#########################
load("/PHENOTYPE/5_lifestyle_v11.RDATA")
source("/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")

# write.table(as.data.frame(PHENO.ANY_SN$sjlid), "/PHENOTYPE/SJLIDs_in_AF_project.txt", col.names = T, quote = F, row.names = F, sep = "\t")

ALL.LIFESTYLE <- edit_lifestyle(ALL.LIFESTYLE)


#########################
## Subsequent Neoplasm ##
#########################

library(haven)
library(benchmarkme)
library(dplyr)
library(plyr)
library(data.table)
library (birk)
library(gtools)
library(stringr)
# library(tidyverse)
library(lubridate)
# benchmarkme::get_ram()
library(survival)

subneo <- read_sas('Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Event Data/subneoplasms.sas7bdat')
common.cols <- colnames(PHENO.ANY_SN)[colnames(PHENO.ANY_SN) %in% colnames(subneo)]
common.cols <- common.cols[-1]
subneo2 <- subneo[,!colnames(subneo) %in% common.cols]
merged_data <- merge(PHENO.ANY_SN, subneo2, by = "sjlid", all.x = TRUE)
merged_data <- merged_data[,!grepl("MRN", colnames(merged_data))]
dim(merged_data)
## Send this to Qi
saveRDS(merged_data, file = "/PHENOTYPE/SJLIFE_complete_data.rds")



head(subneo)
table(subneo$diaggrp)
dim(subneo)
# 1731 9

subneo <- subneo[subneo$sjlid %in% PHENO.ANY_SN$sjlid ,]
dim(subneo)
# 1717
# add diagnosis date 
subneo$diagdt <-  PHENO.ANY_SN$diagdt [match(subneo$sjlid , PHENO.ANY_SN$sjlid)]
subneo$agedx <-  PHENO.ANY_SN$agedx [match(subneo$sjlid , PHENO.ANY_SN$sjlid)]
# add DOB
subneo$DOB <- PHENO.ANY_SN$dob[match(subneo$sjlid, PHENO.ANY_SN$sjlid)]

subneo$gradedt <- as.Date(subneo$gradedt, "%m/%d/%Y") ## **
subneo$AGE.ANY_SN <- time_length(interval(as.Date(subneo$DOB), as.Date(subneo$gradedt)), "years")

subneo$AGE.ANY_SN.after.childhood.cancer <- time_length(interval(as.Date(subneo$diagdt), as.Date(subneo$gradedt)), "years")
subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <- subneo$AGE.ANY_SN - subneo$agedx


########################################
# How many SNs after 5 years
subneo.after5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx > 5,]
length(unique(subneo.after5$sjlid))
# 612

subneo.within5 <- subneo[subneo$AGE.ANY_SN.after.childhood.cancer.from.agedx <= 5,]
sum(!duplicated(subneo.within5$sjlid))
# 22
#############
## Any SNs ##
#############
# Get SNs for the first time and Age at First SN. Note: I am only keeping first event SN 5 years post diagnosis
# For this, I will first sort the table by date
library(data.table)

## Remove SNs as cases that are within 5 years of primary diagnosis
subneo <- subneo[!subneo$sjlid %in% subneo.within5$sjlid,]
ANY_SNs <- setDT(subneo)[,.SD[which.min(gradedt)],by=sjlid][order(gradedt, decreasing = FALSE)]


PHENO.ANY_SN$gradedt <- ANY_SNs$gradedt[match(PHENO.ANY_SN$sjlid, ANY_SNs$sjlid)]
PHENO.ANY_SN$AGE.ANY_SN <- ANY_SNs$AGE.ANY_SN [match(PHENO.ANY_SN$sjlid, ANY_SNs$sjlid)]

# ## Remove SNs if younger than 18 **
# if(sum(PHENO.ANY_SN$AGE.ANY_SN < 18, na.rm = T) > 0){
# PHENO.ANY_SN <- PHENO.ANY_SN[-which(PHENO.ANY_SN$AGE.ANY_SN < 18),]
# }

## remove within 5 years of diagnosis
sum(PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid)
# 22
PHENO.ANY_SN <- PHENO.ANY_SN[!PHENO.ANY_SN$sjlid %in% subneo.within5$sjlid,]


PHENO.ANY_SN$ANY_SNs <- factor(ifelse(!PHENO.ANY_SN$sjlid %in% ANY_SNs$sjlid, 0, 1))

table(PHENO.ANY_SN$ANY_SNs)
# 0    1 
# 3774  605 


#################
## missingness ##
#################
PHENO.ANY_SN$any_tx_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category", "epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_tx_missing  <- factor(ifelse(PHENO.ANY_SN$any_tx_missing == FALSE, "No", "Yes"))

PHENO.ANY_SN$any_chemo_missing <- apply(PHENO.ANY_SN[c("epitxn_dose_5.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_chemo_missing  <- factor(ifelse(PHENO.ANY_SN$any_chemo_missing == FALSE, "No", "Yes"))

PHENO.ANY_SN$any_rt_missing <- apply(PHENO.ANY_SN[c("maxsegrtdose.category", "maxabdrtdose.category", "maxchestrtdose.category")], 1, function(x) any("Unknown" %in% x))
PHENO.ANY_SN$any_rt_missing  <- factor(ifelse(PHENO.ANY_SN$any_rt_missing == FALSE, "No", "Yes"))

PHENO.ANY_SN$any_tx_missing <- factor(PHENO.ANY_SN$any_tx_missing, levels = c("No", "Yes"))
PHENO.ANY_SN$any_chemo_missing <- factor(PHENO.ANY_SN$any_chemo_missing, levels = c("No", "Yes"))
PHENO.ANY_SN$any_rt_missing <- factor(PHENO.ANY_SN$any_rt_missing, levels = c("No", "Yes"))

#########################
## Extract Ethnicities ##
#########################
## Add admixture ethnicity 
ethnicity.admixture <- read.table("/admixture/merged.ancestry.file.txt", header = T)
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, ethnicity.admixture[match(PHENO.ANY_SN$sjlid, ethnicity.admixture$INDIVIDUAL), c("EUR", "EAS", "AFR")])


############################################################
## Drop Unknown level from the lifestyle factor variables ##
############################################################
## Missing tx
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxsegrtdose.category <- droplevels(PHENO.ANY_SN$maxsegrtdose.category)

PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxneckrtdose.category <- droplevels(PHENO.ANY_SN$maxneckrtdose.category)

PHENO.ANY_SN$maxabdrtdose.category[PHENO.ANY_SN$maxabdrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxabdrtdose.category <- droplevels(PHENO.ANY_SN$maxabdrtdose.category)

PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxchestrtdose.category <- droplevels(PHENO.ANY_SN$maxchestrtdose.category)

PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == "Unknown"] <- "None"
PHENO.ANY_SN$maxpelvisrtdose.category <- droplevels(PHENO.ANY_SN$maxpelvisrtdose.category)

PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$epitxn_dose_5.category <- droplevels(PHENO.ANY_SN$epitxn_dose_5.category)

PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$anthra_jco_dose_5.category <- droplevels(PHENO.ANY_SN$anthra_jco_dose_5.category)

PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "Unknown"] <- "None"
PHENO.ANY_SN$aa_class_dose_5.category <- droplevels(PHENO.ANY_SN$aa_class_dose_5.category)


################
## Cross tabs ##
################
library(expss)

# Getting counts for non-missing data only; 6 samples do not have admixture ancestry
CROSS_CASES.df <- PHENO.ANY_SN[!is.na(PHENO.ANY_SN$EUR),]

CROSS_CASES.df <- PHENO.ANY_SN

CROSS_CASES.df <- CROSS_CASES.df[,c("ANY_SNs", "maxsegrtdose.category", "maxchestrtdose.category", "maxabdrtdose.category", "epitxn_dose_5.category")]

CROSS_CASES.df <- apply_labels(CROSS_CASES.df, ANY_SNs = "ANY_SNs", 
                               maxsegrtdose.category = "maxsegrtdose.category", maxchestrtdose.category = "maxchestrtdose.category", 
                               maxabdrtdose.category = "maxabdrtdose.category", epitxn_dose_5.category = "epitxn_dose_5.category")

as.data.frame(t(CROSS_CASES.df %>%
                  cross_cases(ANY_SNs, list(maxsegrtdose.category, maxchestrtdose.category, maxabdrtdose.category, epitxn_dose_5.category))))


cc <- as.data.frame(t(CROSS_CASES.df %>%
                        cross_cases(ANY_SNs, list(maxsegrtdose.category, maxchestrtdose.category, maxabdrtdose.category, epitxn_dose_5.category))))

rownames(cc) <- NULL 
cc


# Create a cross-tabulation table between maxsegrtdose and maxchedtrtdose for cases
cases_table <- table(Max_SegmentedRT_Dose = PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$ANY_SNs == 1],
                     Max_ChestRT_Dose = PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$ANY_SNs == 1])

# Create a cross-tabulation table between maxsegrtdose and maxchedtrtdose for controls
control_table <- table(Max_SegmentedRT_Dose = PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$ANY_SNs == 0],
                     Max_ChestRT_Dose = PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$ANY_SNs == 0])

save.image("analysis_ready.RData")

