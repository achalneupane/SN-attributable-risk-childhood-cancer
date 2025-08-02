
# Load required libraries
library(haven)
library(dplyr)
library(data.table)
library(stringr)
library(lubridate)
library(birk)
library(gtools)
library(benchmarkme)
library(plyr)  # If not used later, consider removing

#########################
## Load Phenotype Data ##
#########################

# Load previous phenotype files
load("/attr_fraction/PHENOTYPE/5_lifestyle_v11_modified_for_HEI_tertiles.RDATA")
pheno_any_sn_old <- ALL.LIFESTYLE
rm(list = setdiff(ls(), "pheno_any_sn_old"))

load("/attr_fraction/PHENOTYPE/3_PRS_scores_categories_v11.RDATA")

#######################################################################
## Recode Physical Activity and Smoking Variables (Based on SAS code) ##
#######################################################################

# Read adult health habits dataset
adult_habits <- read_sas('/Survey Data/adult_healthhabits.sas7bdat')

# Assign participant IDs and remove duplicates
adult_habits <- distinct(adult_habits)
adult_habits$sjlid <- adult_habits$SJLIFEID

# Add Date of Birth from phenotype data
adult_habits$DOB <- pheno_any_sn_old$dob[match(adult_habits$SJLIFEID, pheno_any_sn_old$sjlid)]

# Calculate age at survey; fallback to agesurvey if NA
adult_habits$agesurvey2 <- time_length(interval(as.Date(adult_habits$DOB), as.Date(adult_habits$datecomp)), "years")
adult_habits$agesurvey <- ifelse(is.na(adult_habits$agesurvey2), adult_habits$agesurvey, adult_habits$agesurvey2)

# Filter adult habits to participants present in phenotype data
lifestyle <- adult_habits %>% 
  filter(SJLIFEID %in% pheno_any_sn_old$sjlid) %>%
  mutate(sjlid = SJLIFEID)

##############################
## 1. Physical Activity Code ##
##############################

# vpa10 fixes
lifestyle <- lifestyle %>%
  mutate(
    vpadays = ifelse(vpa10 == 1 & is.na(vpadays), 1, vpadays),
    vpamin = ifelse(vpa10 == 1 & is.na(vpamin), 10, vpamin),
    vpamin = pmin(vpamin, 360)  # Cap at 360 minutes (6 hours)
  )

# Calculate weighted vigorous physical activity (wvpa)
lifestyle <- lifestyle %>%
  mutate(
    wvpa = case_when(
      vpa10 == 1 ~ vpadays * vpamin,
      vpa10 == 2 ~ 0,
      TRUE ~ NA_real_
    )
  )

# Impute wvpa for missing values with conditions
lifestyle <- lifestyle %>%
  mutate(
    wvpa = case_when(
      is.na(wvpa) & (is.na(vpa10) | vpa10 == 2) & nopa == 1 & (is.na(pa20) | pa20 == 0) ~ 0,
      is.na(wvpa) & (is.na(vpa10) | vpa10 == 2) & !is.na(pa20) & pa20 != 0 ~ pa20 * 20,
      is.na(wvpa) & (is.na(vpa10) | vpa10 == 2) & nopa == 2 & (is.na(pa20) | pa20 == 0) ~ 0,
      TRUE ~ wvpa
    )
  )

# mpa10 fixes
lifestyle <- lifestyle %>%
  mutate(
    mpa10 = ifelse(is.na(mpa10) & (!is.na(mpadays) | !is.na(mpamin)), 1, mpa10),
    mpadays = ifelse(mpa10 == 1 & is.na(mpadays), 1, mpadays),
    mpamin = ifelse(mpa10 == 1 & is.na(mpamin), 1, mpamin),
    mpamin = pmin(mpamin, 360)  # Cap at 360 minutes
  )

# Calculate weighted moderate physical activity (wmpa)
lifestyle <- lifestyle %>%
  mutate(
    wmpa = case_when(
      mpa10 == 1 ~ mpadays * mpamin,
      mpa10 == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    wmpa = ifelse(is.na(wmpa) & (is.na(mpa10) | mpa10 == 2) & nopa == 1, 0, wmpa),
    wmpa = ifelse(is.na(wmpa) & !is.na(wvpa), 0, wmpa)
  )

# Calculate combined weighted physical activity with cap
lifestyle <- lifestyle %>%
  mutate(
    mvpawk = wmpa + (wvpa * 2),
    mvpawk = pmin(mvpawk, 2520),  # Cap at 2520 minutes (6 hours per day)
    CDC_PA = ifelse(mvpawk < 150, 0, 1)
  )



## 1. New Physical Activity (age ≥ 18)
MET_iid_dob_18 <- subset(lifestyle, agesurvey >= 18)
MET_iid_dob_18 <- MET_iid_dob_18[!is.na(MET_iid_dob_18$mvpawk), ]
MET_iid_dob_18 <- MET_iid_dob_18[order(MET_iid_dob_18$sjlid, MET_iid_dob_18$agesurvey), ]
MET_iid_dob_18_uniq <- MET_iid_dob_18[!duplicated(MET_iid_dob_18$sjlid), ]
MET_iid_dob_18_uniq$PhysicalActivity_yn <- MET_iid_dob_18_uniq$mvpawk


## 2. Smoking Status - merge variables from adult_habbits into lifestyle
lifestyle$evsm <- adult_habbits$evsm[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$smnow <- adult_habbits$smnow[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$cigmo <- adult_habbits$cigmo[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$cigd <- adult_habbits$cigd[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]
lifestyle$smyr <- adult_habbits$smyr[match(lifestyle$SJLIFEID, adult_habbits$SJLIFEID)]

# Fix missing or inconsistent smoking variables based on evsm
lifestyle$smnow[lifestyle$evsm == 2] <- 2
lifestyle$cigmo[lifestyle$evsm == 2] <- 2
lifestyle$cigd[lifestyle$evsm == 2] <- 0
lifestyle$smyr[lifestyle$evsm == 2] <- 0

lifestyle$evsm[is.na(lifestyle$evsm) & lifestyle$smnow == 1] <- 1
lifestyle$smnow[lifestyle$evsm == 1 & is.na(lifestyle$smnow) & lifestyle$cigmo == 1] <- 1
lifestyle$smnow[lifestyle$evsm == 1 & is.na(lifestyle$smnow) & lifestyle$cigmo == 2] <- 2

# Create smoking status categories: 1=Former, 2=Current, 3=Never
lifestyle$smkStat <- NA_integer_
lifestyle$smkStat[lifestyle$evsm == 1 & lifestyle$smnow == 2] <- 1
lifestyle$smkStat[lifestyle$evsm == 1 & lifestyle$smnow == 1] <- 2
lifestyle$smkStat[lifestyle$evsm == 2] <- 3

table(lifestyle$smkStat)


## Smoking subset for adults (age ≥ 18)
smk_iid_dob_18 <- subset(lifestyle, agesurvey >= 18)
smk_iid_dob_18 <- smk_iid_dob_18[!is.na(smk_iid_dob_18$cigd), ]
smk_iid_dob_18 <- smk_iid_dob_18[order(smk_iid_dob_18$sjlid, smk_iid_dob_18$agesurvey), ]
smk_iid_dob_18_uniq <- smk_iid_dob_18[!duplicated(smk_iid_dob_18$sjlid), ]

smk_iid_dob_18_uniq$smoker_former_or_never_yn <- smk_iid_dob_18_uniq$cigd
table(smk_iid_dob_18_uniq$cigd)


## 3. Drinking status (age ≥ 18)
drk_iid_dob_18 <- subset(lifestyle, agesurvey >= 18)
drk_iid_dob_18$RiskyHeavyDrink_yn <- drk_iid_dob_18$frdrk
drk_iid_dob_18 <- drk_iid_dob_18[!is.na(drk_iid_dob_18$RiskyHeavyDrink_yn), ]
drk_iid_dob_18 <- drk_iid_dob_18[order(drk_iid_dob_18$sjlid, drk_iid_dob_18$agesurvey), ]
drk_iid_dob_18_uniq <- drk_iid_dob_18[!duplicated(drk_iid_dob_18$sjlid), ]


## 4. Adult BMI (using Yadav’s data with more samples)
library(sas7bdat)

bmi <- read_sas('/Clinical Data/function_combo_basic.sas7bdat')
iid <- read.table('/attr_fraction/4401_attributable_fraction_ids.txt', header = FALSE)

bmi_iid <- subset(bmi, sjlid %in% iid$V1, select = c('sjlid', 'assmntdate', 'BMIadj'))
demo <- read_sas('/Clinical Data/demographics.sas7bdat')
dob <- demo[c('sjlid', 'dob')]

bmi_iid_dob <- merge(bmi_iid, dob, by = 'sjlid')
bmi_iid_dob$agebmi <- as.numeric(difftime(as.Date(bmi_iid_dob$assmntdate), as.Date(bmi_iid_dob$dob), units = "days")) / 365.25

bmi_iid_dob_18 <- subset(bmi_iid_dob, agebmi >= 18)
bmi_iid_dob_18 <- bmi_iid_dob_18[order(bmi_iid_dob_18$sjlid, bmi_iid_dob_18$agebmi), ]
bmi_iid_dob_18_uniq <- bmi_iid_dob_18[!duplicated(bmi_iid_dob_18$sjlid), ]

bmi_iid_dob_18_uniq$Not_obese_yn <- as.numeric(bmi_iid_dob_18_uniq$BMIadj < 30)


## Part 2. Adult habits from Siddhant
adlthabits <- read.delim("/attr_fraction/PHENOTYPE/adlthabits.txt", header = TRUE, sep = "\t")
adlthabits$sjlid <- adlthabits$SJLIFEID

# Fix date format (MM/DD/YYYY to YYYY-MM-DD)
adlthabits$datecomp <- paste(
  sapply(strsplit(adlthabits$datecomp, "/"), `[`, 3),
  sapply(strsplit(adlthabits$datecomp, "/"), `[`, 1),
  sapply(strsplit(adlthabits$datecomp, "/"), `[`, 2),
  sep = "-"
)

adlthabits <- distinct(adlthabits)

# Get DOB from PHENO.ANY_SN and calculate agesurvey
adlthabits$DOB <- PHENO.ANY_SN$dob[match(adlthabits$SJLIFEID, PHENO.ANY_SN$sjlid)]
library(lubridate)
adlthabits$agesurvey2 <- time_length(interval(as.Date(adlthabits$DOB), as.Date(adlthabits$datecomp)), "years")
adlthabits$agesurvey2[is.na(adlthabits$agesurvey2)] <- adlthabits$agesurvey[is.na(adlthabits$agesurvey2)]
adlthabits$agesurvey <- adlthabits$agesurvey2

adlthabits <- adlthabits[adlthabits$sjlid %in% PHENO.ANY_SN$sjlid, ]


## Smoking status from Siddhant data (age ≥ 18)
smk_iid_dob_18.2 <- subset(adlthabits, agesurvey >= 18)
smk_iid_dob_18.2 <- smk_iid_dob_18.2[!is.na(smk_iid_dob_18.2$smoker), ]
smk_iid_dob_18.2 <- smk_iid_dob_18.2[order(smk_iid_dob_18.2$sjlid, smk_iid_dob_18.2$agesurvey), ]
smk_iid_dob_18_uniq.2 <- smk_iid_dob_18.2[!duplicated(smk_iid_dob_18.2$sjlid), ]

# Recode smoker status: former/never = 1, current = 0 (binary)
smk_iid_dob_18_uniq.2$smoker_former_or_never_yn <- as.numeric(smk_iid_dob_18_uniq.2$smoker != 2)
smk_iid_dob_18_uniq.2$smoker_ever_yn <- ifelse(smk_iid_dob_18_uniq.2$smoker != 3, "Yes", "No")

table(smk_iid_dob_18_uniq.2$smoker_ever_yn)
table(smk_iid_dob_18_uniq.2$smoker)
table(smk_iid_dob_18_uniq$smkStat)


## Drinking status from Siddhant data (age ≥ 18)
# Replace zeros in binge, heavy, risky drink columns with 0 to ensure numeric consistency
cols_drink <- grep("bingedrink|heavydrink|riskydrink", colnames(adlthabits), value = TRUE)
adlthabits[cols_drink] <- lapply(adlthabits[cols_drink], function(x) ifelse(x == 0, 0, x))

drk_iid_dob_18.2 <- subset(adlthabits, agesurvey >= 18)
# Define NOT risky heavy drink = 1 if both heavy and risky drinking = 0
drk_iid_dob_18.2$Not_risky_heavy_drink <- as.numeric(
  (drk_iid_dob_18.2$heavydrink == 0 & drk_iid_dob_18.2$riskydrink == 0)
)

table(drk_iid_dob_18.2$Not_risky_heavy_drink)


## Dietary habits from Siddhant data (age ≥ 18)
diet_iid_dob_18.2 <- subset(adlthabits, agesurvey >= 18)
diet_iid_dob_18.2 <- diet_iid_dob_18.2[order(diet_iid_dob_18.2$sjlid, diet_iid_dob_18.2$agesurvey), ]
diet_iid_dob_18_uniq.2 <- diet_iid_dob_18.2[!duplicated(diet_iid_dob_18.2$sjlid), ]


##########################################
## Section 5: Extract and Clean Diet Data
##########################################

# Load earliest diet data after age 18
adultdiet <- read.table("/attr_fraction/PHENOTYPE/adultbmi.txt", sep = "\t", header = TRUE)
adultdiet <- adultdiet[adultdiet$DateVisitStart != "", ]

# Format DateVisitStart to YYYY-MM-DD
adultdiet$DateVisitStart <- with(
  strsplit(adultdiet$DateVisitStart, "\\/"),
  paste(sapply(., `[`, 3), sapply(., `[`, 1), sapply(., `[`, 2), sep = "-")
)

# Load additional dietary variables from original SAS dataset
original.adultdiet <- read_sas("/Clinical Data/ffq_grams.sas7bdat")
original.adultdiet$DateVisitStart <- gsub("-0", "-", original.adultdiet$DateVisitStart)

# Create a KEY variable for matching
original.adultdiet$KEY <- paste(original.adultdiet$sjlid, original.adultdiet$DateVisitStart, sep = ":")
adultdiet$KEY <- paste(adultdiet$sjlid, adultdiet$DateVisitStart, sep = ":")

# Match and align original variables
matched_keys <- adultdiet$KEY %in% original.adultdiet$KEY
original.adultdiet <- original.adultdiet[original.adultdiet$KEY %in% adultdiet$KEY, ]
original.adultdiet <- original.adultdiet[match(adultdiet$KEY, original.adultdiet$KEY), ]

# Add missing variables from original data
add_vars <- c("NUTSFREQ", "SOFTDRINKSFREQ", "MIXEDBEEFPORKFREQ", "NOTFRIEDFISHFREQ",
              "HOTDOGFREQ", "BACONFREQ", "SAUSAGEFREQ", "G_NWHL")
adultdiet[add_vars] <- lapply(original.adultdiet[add_vars], as.numeric)

# Replace AGE with calculated value from DOB
adultdiet$AGE <- as.numeric(adultdiet$AGE)
adultdiet$DOB <- PHENO.ANY_SN$dob[match(adultdiet$sjlid, PHENO.ANY_SN$sjlid)]
adultdiet$AGE_at_Visit <- time_length(interval(as.Date(adultdiet$DOB), as.Date(adultdiet$DateVisitStart)), "years")

# Keep unique participants in PHENO.ANY_SN
samples.sjlife <- intersect(unique(adultdiet$sjlid), PHENO.ANY_SN$sjlid)

##################
## Define DIET ##
##################

DIET <- adultdiet

# Binary dietary indicators
DIET$FRUITSRV_yn               <- as.numeric(DIET$FRUITSRV >= 3)
DIET$NUTSFREQ_yn               <- as.numeric(DIET$NUTSFREQ >= 5)
DIET$VEGSRV_yn                 <- as.numeric(DIET$VEGSRV >= 3)
DIET$WGRAINS_yn                <- as.numeric(DIET$WGRAINS >= 3)
DIET$NOTFRIEDFISHFREQ_yn       <- as.numeric(DIET$NOTFRIEDFISHFREQ >= 6)
DIET$DAIRYSRV_yn               <- as.numeric(DIET$DAIRYSRV >= 2.5)
DIET$MIXEDBEEFPORKFREQ_yn      <- as.numeric(DIET$MIXEDBEEFPORKFREQ <= 5)
DIET$DT_TFAT_cohort_median_yn  <- as.numeric(DIET$DT_TFAT <= median(DIET$DT_TFAT, na.rm = TRUE))
DIET$SOFTDRINKSFREQ_yn         <- as.numeric(DIET$SOFTDRINKSFREQ <= 5)
DIET$DT_SODI_yn                <- as.numeric(DIET$DT_SODI <= 2000)

# Healthy diet: meets ≥5 of the 10 criteria
diet_cols <- c("FRUITSRV_yn", "NUTSFREQ_yn", "VEGSRV_yn", "WGRAINS_yn", "NOTFRIEDFISHFREQ_yn", 
               "DAIRYSRV_yn", "MIXEDBEEFPORKFREQ_yn", "DT_TFAT_cohort_median_yn", 
               "SOFTDRINKSFREQ_yn", "DT_SODI_yn")
DIET$HEALTHY_Diet_yn <- as.numeric(rowSums(DIET[, diet_cols], na.rm = TRUE) >= 5)

# Use calculated AGE if AGE_at_Visit is NA
DIET$AGE_at_Visit[is.na(DIET$AGE_at_Visit)] <- DIET$AGE[is.na(DIET$AGE_at_Visit)]
DIET$agesurvey <- as.numeric(DIET$AGE_at_Visit)

###################################################
## Select First Diet Record After Age 18 per ID ##
###################################################

filter_first_record <- function(df, score_var) {
  # Filter to age 18+ and non-missing score_var
  df_filtered <- df[df$agesurvey >= 18 & !is.na(df[[score_var]]), ]
  # Order by ID and survey age
  df_filtered <- df_filtered[order(df_filtered$sjlid, df_filtered$agesurvey), ]
  # Keep only first record per ID
  df_filtered[!duplicated(df_filtered$sjlid), ]
}

# Apply function to diet scores
diet_iid_dob_18_uniq    <- filter_first_record(DIET, "HEALTHY_Diet_yn")
HEI2005_iid_dob_18_uniq <- filter_first_record(DIET, "HEI2005_TOTAL_SCORE")
HEI2010_iid_dob_18_uniq <- filter_first_record(DIET, "HEI2010_TOTAL_SCORE")
HEI2015_iid_dob_18_uniq <- filter_first_record(DIET, "HEI2015_TOTAL_SCORE")

# Final ID list from phenotype data
sjlid <- PHENO.ANY_SN$sjlid

# Helper function to safely match and extract columns
match_column <- function(ids, df, col_name) {
  df[[col_name]][match(ids, df$sjlid %||% df$SJLIFEID)]
}

# Initialize ALL.LIFESTYLE with SJLIFEID and PhysicalActivity data
ALL.LIFESTYLE <- data.frame(
  SJLIFEID = sjlid,
  PhysicalActivity_yn = match_column(sjlid, MET_iid_dob_18_uniq, "PhysicalActivity_yn"),
  PhysicalActivity_yn_agesurvey = match_column(sjlid, MET_iid_dob_18_uniq, "agesurvey")
)

# Add smoker variables
ALL.LIFESTYLE$smoker_former_or_never_yn <- match_column(sjlid, smk_iid_dob_18_uniq, "smoker_former_or_never_yn")
ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey <- match_column(sjlid, smk_iid_dob_18_uniq, "agesurvey")

# Duplicate smoker variable under a different name
ALL.LIFESTYLE$Smoker_ever_yn <- ALL.LIFESTYLE$smoker_former_or_never_yn
ALL.LIFESTYLE$Smoker_ever_yn_agesurvey <- ALL.LIFESTYLE$smoker_former_or_never_yn_agesurvey

# Add drinking variables
ALL.LIFESTYLE$RiskyHeavyDrink_yn <- match_column(sjlid, drk_iid_dob_18_uniq, "RiskyHeavyDrink_yn")
ALL.LIFESTYLE$RiskyHeavyDrink_yn_agesurvey <- match_column(sjlid, drk_iid_dob_18_uniq, "agesurvey")

# Add obesity variables
ALL.LIFESTYLE$Obese_yn <- match_column(sjlid, bmi_iid_dob_18_uniq, "BMIadj")
ALL.LIFESTYLE$Obese_yn_agesurvey <- match_column(sjlid, bmi_iid_dob_18_uniq, "agebmi")

# Optionally, print frequency tables for quick checks
print(table(ALL.LIFESTYLE$PhysicalActivity_yn))
print(table(ALL.LIFESTYLE$smoker_former_or_never_yn))
print(table(ALL.LIFESTYLE$RiskyHeavyDrink_yn))
print(table(ALL.LIFESTYLE$Obese_yn))

# Assign diet variables with matching and factor conversion for HEALTHY_Diet_yn
ALL.LIFESTYLE$HEALTHY_Diet_yn <- diet_iid_dob_18_uniq$HEALTHY_Diet_yn[match(ALL.LIFESTYLE$SJLIFEID, diet_iid_dob_18_uniq$sjlid)]
ALL.LIFESTYLE$HEALTHY_Diet_yn_agesurvey <- diet_iid_dob_18_uniq$agesurvey[match(ALL.LIFESTYLE$SJLIFEID, diet_iid_dob_18_uniq$sjlid)]

ALL.LIFESTYLE$HEALTHY_Diet_yn[is.na(ALL.LIFESTYLE$HEALTHY_Diet_yn)] <- "Unknown"
ALL.LIFESTYLE$HEALTHY_Diet_yn <- factor(ALL.LIFESTYLE$HEALTHY_Diet_yn, levels = c(1, 0, "Unknown"))
table(ALL.LIFESTYLE$HEALTHY_Diet_yn)
# 1       0 Unknown 
# 306    2850    1245

# Assign HEI diet scores and ages
for (score in c("HEI2005_TOTAL_SCORE", "HEI2010_TOTAL_SCORE", "HEI2015_TOTAL_SCORE")) {
  uniq_obj <- get(paste0(score, "_iid_dob_18_uniq"))
  ALL.LIFESTYLE[[score]] <- uniq_obj[[score]][match(ALL.LIFESTYLE$SJLIFEID, uniq_obj$sjlid)]
  ALL.LIFESTYLE[[paste0(score, "_agesurvey")]] <- uniq_obj$agesurvey[match(ALL.LIFESTYLE$SJLIFEID, uniq_obj$sjlid)]
}

#########################
## Create HEI tertiles ##
#########################
keep.reference <- PHENO.ANY_SN.old$SJLIFEID[PHENO.ANY_SN.old$HEI2015_TOTAL_SCORE.lt60.category == "No"]

# Prepare tertile categories, set to NA for references
ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category <- ALL.LIFESTYLE$HEI2015_TOTAL_SCORE
ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category[ALL.LIFESTYLE$SJLIFEID %in% keep.reference] <- NA

# Calculate tertile breaks excluding zero and NA
tertiles <- unname(quantile(ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category[ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category != 0], probs = c(1/3, 2/3, 1), na.rm = TRUE))

ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category <- cut(
  ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category,
  breaks = c(0, tertiles),
  labels = c("3rd", "2nd", "1st"),
  include.lowest = TRUE
)

ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category[ALL.LIFESTYLE$SJLIFEID %in% keep.reference] <- "No"
ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category[is.na(ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category)] <- "Unknown"

ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category <- factor(
  ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category,
  levels = c("No", "3rd", "2nd", "1st", "Unknown")
)

table(ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category)

## Physical activity tertiles
ALL.LIFESTYLE$PhysicalActivity_yn <- as.numeric(ALL.LIFESTYLE$PhysicalActivity_yn)

keep.reference <- PHENO.ANY_SN.old$SJLIFEID[PHENO.ANY_SN.old$PhysicalActivity_yn == 1]
ALL.LIFESTYLE$PhysicalActivity_yn[ALL.LIFESTYLE$SJLIFEID %in% keep.reference] <- NA

tertiles <- unname(quantile(ALL.LIFESTYLE$PhysicalActivity_yn[ALL.LIFESTYLE$PhysicalActivity_yn != 0], probs = c(1/3, 2/3, 1), na.rm = TRUE))

ALL.LIFESTYLE$PhysicalActivity_yn <- cut(
  ALL.LIFESTYLE$PhysicalActivity_yn,
  breaks = c(0, tertiles),
  labels = c("3rd", "2nd", "1st"),
  include.lowest = TRUE
)

ALL.LIFESTYLE$PhysicalActivity_yn[ALL.LIFESTYLE$SJLIFEID %in% keep.reference] <- "Yes"
ALL.LIFESTYLE$PhysicalActivity_yn[is.na(ALL.LIFESTYLE$PhysicalActivity_yn)] <- "Unknown"

ALL.LIFESTYLE$PhysicalActivity_yn <- factor(
  ALL.LIFESTYLE$PhysicalActivity_yn,
  levels = c("Yes", "3rd", "2nd", "1st", "Unknown")
)

table(ALL.LIFESTYLE$PhysicalActivity_yn)

## Smoking tertiles
table(PHENO.ANY_SN.old$Smoker_ever_yn)
# No     Yes Unknown 
# 2425    1129       0 

ALL.LIFESTYLE$Smoker_ever_yn <- as.numeric(ALL.LIFESTYLE$Smoker_ever_yn)
keep.reference <- PHENO.ANY_SN.old$SJLIFEID[PHENO.ANY_SN.old$Smoker_ever_yn == "No"]

tertiles <- unname(quantile(ALL.LIFESTYLE$Smoker_ever_yn[ALL.LIFESTYLE$Smoker_ever_yn != 0], probs = c(1/3, 2/3, 1), na.rm = TRUE))

ALL.LIFESTYLE$Smoker_ever_yn <- cut(
  ALL.LIFESTYLE$Smoker_ever_yn,
  breaks = c(0, tertiles),
  labels = c("1st", "2nd", "3rd"),
  include.lowest = TRUE
)

ALL.LIFESTYLE$Smoker_ever_yn[ALL.LIFESTYLE$SJLIFEID %in% keep.reference] <- "No"
ALL.LIFESTYLE$Smoker_ever_yn[is.na(ALL.LIFESTYLE$Smoker_ever_yn)] <- "Unknown"

ALL.LIFESTYLE$Smoker_ever_yn <- factor(
  ALL.LIFESTYLE$Smoker_ever_yn,
  levels = c("No", "1st", "2nd", "3rd", "Unknown")
)

table(ALL.LIFESTYLE$Smoker_ever_yn)

## Drinking tertiles
keep.reference <- PHENO.ANY_SN.old$SJLIFEID[PHENO.ANY_SN.old$NOT_RiskyHeavyDrink_yn == 1]

ALL.LIFESTYLE$RiskyHeavyDrink_yn <- as.numeric(ALL.LIFESTYLE$RiskyHeavyDrink_yn)
ALL.LIFESTYLE$RiskyHeavyDrink_yn[ALL.LIFESTYLE$SJLIFEID %in% keep.reference] <- NA

tertiles <- unname(quantile(ALL.LIFESTYLE$RiskyHeavyDrink_yn[ALL.LIFESTYLE$RiskyHeavyDrink_yn != 0], probs = c(1/3, 2/3, 1), na.rm = TRUE))

ALL.LIFESTYLE$RiskyHeavyDrink_yn <- cut(
  ALL.LIFESTYLE$RiskyHeavyDrink_yn,
  breaks = c(0, tertiles),
  labels = c("1st", "2nd", "3rd"),
  include.lowest = TRUE
)

ALL.LIFESTYLE$RiskyHeavyDrink_yn[is.na(ALL.LIFESTYLE$RiskyHeavyDrink_yn)] <- "Unknown"
ALL.LIFESTYLE$RiskyHeavyDrink_yn[ALL.LIFESTYLE$SJLIFEID %in% keep.reference] <- "No"

ALL.LIFESTYLE$RiskyHeavyDrink_yn <- factor(
  ALL.LIFESTYLE$RiskyHeavyDrink_yn,
  levels = c("No", "1st", "2nd", "3rd", "Unknown")
)

table(ALL.LIFESTYLE$RiskyHeavyDrink_yn)
# --- Obesity variable recoding ---
# Recode BMI categories according to CDC adult BMI categories:
# Underweight (<18.5), Normal weight (18.5-24.9), Overweight (25-29.9), Obese (>=30)
ALL.LIFESTYLE$Obese_yn <- cut(
  ALL.LIFESTYLE$Obese_yn,
  breaks = c(-Inf, 18.5, 24.9, 29.9, Inf),
  labels = c("Underweight", "No", "Overweight", "Obese"),
  right = FALSE
)

# Reorder factor levels so "No" (normal weight) is the reference level
ALL.LIFESTYLE$Obese_yn <- factor(
  ALL.LIFESTYLE$Obese_yn,
  levels = c("No", "Underweight", "Overweight", "Obese", "Unknown")
)

# --- HEI tertiles ---
# Preserve HEI tertiles from the original tertile category
ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.lt60.category <- ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category

# Optional: Create a summary dataframe for plotting or checking HEI scores
hei_summary <- data.frame(
  SJLIFEID = ALL.LIFESTYLE$SJLIFEID,
  HEI_Tertile = ALL.LIFESTYLE$HEI2015_TOTAL_SCORE.tertile.category,
  HEI_Score = ALL.LIFESTYLE$HEI2015_TOTAL_SCORE
)

# --- Load and verify clinical data from 2020 freeze ---
chemo_2020 <- read_sas('/Clinical Data/chemosum_dose.sas7bdat')
radiation_2020 <- read_sas('/Clinical Data/radiation_dosimetry.sas7bdat')
demographics <- read_sas('/Clinical Data/demographics.sas7bdat')
diagnosis <- read_sas('/Clinical Data/diagnosis.sas7bdat')

# --- Merge key variables to PHENO.ANY_SN by SJLIFE ID ---
PHENO.ANY_SN$gender <- factor(demographics$gender[match(PHENO.ANY_SN$sjlid, demographics$sjlid)], levels = c("Male", "Female"))
PHENO.ANY_SN$dob <- demographics$dob[match(PHENO.ANY_SN$sjlid, demographics$sjlid)]
PHENO.ANY_SN$agedx <- diagnosis$agedx[match(PHENO.ANY_SN$sjlid, diagnosis$sjlid)]
PHENO.ANY_SN$diagdt <- diagnosis$diagdt[match(PHENO.ANY_SN$sjlid, diagnosis$sjlid)]

# --- Create age at diagnosis categories ---
PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- cut(
  PHENO.ANY_SN$agedx,
  breaks = c(-Inf, 5, 10, 15, Inf),
  labels = c("0-4", "5-9", "10-14", ">=15"),
  right = FALSE
)
PHENO.ANY_SN$AGE_AT_DIAGNOSIS <- factor(
  PHENO.ANY_SN$AGE_AT_DIAGNOSIS,
  levels = c("0-4", "5-9", "10-14", ">=15")
)


# Load age at last contact and match by sjlid
lstcondt <- read_sas('/Tracking Data/lstcondt.sas7bdat')
PHENO.ANY_SN$agelstcontact <- lstcondt$agelstcontact[match(PHENO.ANY_SN$sjlid, lstcondt$sjlid)]

# Remove unnecessary columns
drop_cols <- c(
  "aa_class_dose_any_yn", "aa_class_dose_any.category", "aa_class_dose_5_yn",
  "cisplat_dose_any_yn", "cisplat_dose_any.category", "cisplateq_dose_5_yn", "cisplateq_dose_5.category",
  "aa_hvymtl_dose_any_yn", "aa_hvymtl_dose_any.category", "aa_hvymtl_dose_5_yn", "aa_hvymtl_dose_5.category",
  "anthra_jco_dose_any_yn", "anthra_jco_dose_any.category", "anthra_jco_dose_5_yn",
  "epitxn_dose_any_yn", "epitxn_dose_any.category", "epitxn_dose_5_yn",
  "aa_class_dose_any", "aa_hvymtl_dose_5", "aa_hvymtl_dose_any",
  "cisplat_dose_5", "cisplat_dose_any", "cisplateq_dose_5", "cisplateq_dose_any",
  "epitxn_dose_any", "anthra_cog_dose_5", "anthra_cog_dose_any",
  "anthra_jco_dose_any", "carbo_dose_any", "carbo_dose_5"
)
PHENO.ANY_SN <- PHENO.ANY_SN[, !colnames(PHENO.ANY_SN) %in% drop_cols]

# Helper function to match dose, cut into categories, and handle unknowns
create_dose_category <- function(df, dose_var, match_df, sjlid_var = "sjlid", 
                                 breaks, labels) {
  doses <- match_df[[dose_var]][match(df[[sjlid_var]], match_df[[sjlid_var]])]
  cat_var <- cut(doses, breaks = breaks, labels = labels, include.lowest = TRUE)
  levels(cat_var) <- c(levels(cat_var), "Unknown")
  cat_var[is.na(doses)] <- "Unknown"
  list(dose = doses, category = cat_var)
}

# Radiation doses and categories
radiation_vars <- list(
  maxsegrtdose = list(breaks = c(0, 200, 1799, 2999, Inf), labels = c("None", ">0-<18", ">=18-<30", ">=30")),
  maxchestrtdose = list(breaks = c(0, 200, 1999, Inf), labels = c("None", ">0-<20", ">=20")),
  maxneckrtdose = list(breaks = c(0, 200, 1099, 1999, 2999, Inf), labels = c("None", ">0-<11", ">=11-<20", ">=20-<30", ">=30")),
  maxabdrtdose = list(breaks = c(0, 200, 2999, Inf), labels = c("None", ">0-<30", ">=30")),
  maxpelvisrtdose = list(breaks = c(0, 200, 1999, Inf), labels = c("None", ">0-<20", ">=20"))
)

for (var in names(radiation_vars)) {
  result <- create_dose_category(
    PHENO.ANY_SN, var, radiation.2020,
    breaks = radiation_vars[[var]]$breaks,
    labels = radiation_vars[[var]]$labels
  )
  PHENO.ANY_SN[[var]] <- result$dose
  PHENO.ANY_SN[[paste0(var, ".category")]] <- result$category
}

# Helper function for chemo doses with tertile categories
create_tertile_category <- function(df, dose_var, match_df, sjlid_var = "sjlid") {
  doses <- match_df[[dose_var]][match(df[[sjlid_var]], match_df[[sjlid_var]])]
  tertiles <- unname(quantile(doses[doses != 0], probs = c(1/3, 2/3, 1), na.rm = TRUE))
  breaks <- c(-Inf, 0, tertiles)
  labels <- c("None", "1st", "2nd", "3rd")
  cat_var <- cut(doses, breaks = breaks, labels = labels, include.lowest = TRUE)
  levels(cat_var) <- c(levels(cat_var), "Unknown")
  cat_var[is.na(cat_var)] <- "Unknown"
  list(dose = doses, category = cat_var)
}

# Chemo dose variables and categories
chemo_vars <- c(
  aa_class_dose_5 = "alkylating_dose_5",
  anthra_jco_dose_5 = "anthracyclines_dose_5",
  epitxn_dose_5 = "epipodophyllotoxins_dose_5"
)

for (pheno_var in names(chemo_vars)) {
  chemo_var <- chemo_vars[[pheno_var]]
  result <- create_tertile_category(PHENO.ANY_SN, chemo_var, chemo.2020)
  PHENO.ANY_SN[[pheno_var]] <- result$dose
  PHENO.ANY_SN[[paste0(pheno_var, ".category")]] <- result$category
}

# Save workspace
save.image("/attr_fraction/PHENOTYPE/lancetOncology/5_lifestyle_v11_modified_for_HEI_tertiles.RDATA")


# load("/attr_fraction/PHENOTYPE/5_lifestyle_v11_modified_for_HEI_tertiles.RDATA")