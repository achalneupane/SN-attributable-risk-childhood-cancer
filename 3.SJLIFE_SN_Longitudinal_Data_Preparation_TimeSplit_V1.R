
################## # Prepare data ## ################# 
#This R code prepares longitudinal clinical data on subsequent neoplasms (SNs)
#in childhood cancer survivors for survival analysis by restructuring the
#dataset into time-to-event intervals. It processes multiple events per
#individual by ordering events, defining start and end times for each risk
#interval beginning 5 years after initial diagnosis, and differentiating event
#and censoring periods. The code handles recurrent events by creating sequential
#intervals between diagnosis, events, and last follow-up, then splits the
#follow-up time into yearly segments to facilitate time-varying hazard modeling.
#Additionally, it generates cubic spline variables for attained age to allow
#flexible adjustment of age effects in downstream survival models. Overall, this
#preprocessing transforms raw event data into a format suitable for advanced
#survival analysis of SN risk over time.

load("analysis_ready.RData")
# PHENO.ANY_SN <- merged_data
data <- PHENO.ANY_SN[!grepl("MRN", colnames(PHENO.ANY_SN))]

## Age at last contact for cases is SN diagnosis data
data$agelstcontact[!is.na(data$AGE.ANY_SN)] <- data$AGE.ANY_SN[!is.na(data$AGE.ANY_SN)]

data$event <- ifelse(!is.na(data$gradedt), 1, 0)

data$first <- ave(data$agelstcontact, data$sjlid, FUN = seq_along)
M <- max(data$first, na.rm = T)  ### maximum number of events
data[data$first==M,]$sjlid  #the id for the person with the maximum number of rows.
data[data$sjlid==data[data$first==M,]$sjlid,]

### Take the last row so we know the maximum number of events per person
event.number <- do.call(rbind, lapply(split(data, data$sjlid), tail, 1))[,c("sjlid","first")]
colnames(event.number) <- c("sjlid","maxE")

alldata <- merge(data,event.number,by.x="sjlid",by.y="sjlid")

alldata$start <- NULL
alldata$end <- NULL
### For those without event, start is agedx and end is Fu date === for SN, analysis starts from 5 years post DX (SNs within 5 years have been removed from the analysis)
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0] + 5
### Achal: Since the analysis start from 5 years post DX, the above line has been revised to: alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0]+5
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0] 

### For the first event, start is agedx and end is first event time
alldata$start[alldata$event==1 & alldata$first==1] <- alldata$agedx[alldata$event==1 & alldata$first==1] + 5
### Achal: also added +5 in the above line
alldata$end[alldata$event==1 & alldata$first==1] <- as.numeric(difftime(alldata$gradedt[alldata$event==1 & alldata$first==1],alldata$dob[alldata$event==1 & alldata$first==1], units = 'days')/365.25)

#### For events that are not the first, segments are from the previous event date to this event date
alldata$previous <- as.Date(c(NA,alldata$gradedt[1:length(alldata$gradedt)-1]),origin = "1970-01-01")
alldata[1:10,]
### if first>1 (i.e, 2 or more events, previous event time remained, others are missing)
alldata$previous[alldata$first==1] <- NA
alldata$start[alldata$first>1] <- as.numeric(difftime(alldata$previous[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)
alldata$end[alldata$first>1] <- as.numeric(difftime(alldata$gradedt[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)


### If one person has only 1 event, we need to have segment from agedx to event date (handled above), and then from event date to end of FU
### If one person has multiple events, we need to add segments from the last event to end of Fu.
row_add <- alldata[alldata$event==1 & alldata$first==alldata$maxE,] ## this includes (1) people with only 1 event (2) The last row for people with multiple events
row_add$start <- as.numeric(difftime(row_add$gradedt, row_add$dob, units = 'days')/365.25)
row_add$end <- row_add$agelstcontact;
#### remember to make these rows event as 0, as we are adding the time after the last event date.
row_add$event <- 0
### If we do first event analysis, we would not need these. I need to add an indicator here, in case I only need to get the first event analysis segments.
row_add$evt1 <- 0;

alldata$evt1 <- 0;
alldata$evt1[alldata$first==1] <- 1
adata <- rbind(alldata,row_add)
#### order by person and start date
adata <- adata[order(adata$sjlid, adata$start, decreasing = FALSE),]
table(adata$event)## Double check event numebr is correct

###any stop time <=start time
adata$end[adata$end<=adata$start] <- adata$end[adata$end<=adata$start] + 1/365
any <- adata[adata$end<=adata$start,] # 
diff=any$start-any$end 
dim(adata)
final <- adata[adata$end>adata$start,]

# final$sjlid[duplicated(final$sjlid)]
# cc.final <- final[, c(1,(ncol(final)-11):ncol(final))]
# cc.data <- data[, c(1,(ncol(data)-11):ncol(data))]

minimum  <-  min(final$start, na.rm = TRUE) - 1
if (minimum < 0) minimum <- 0
maximum  <-  max(final$end, na.rm = TRUE) + 1
SNs_py <-  survSplit(final, cut=seq(minimum, maximum, 1), end="end",start="start",event="event") 
table(SNs_py$event)   ## Double check event numebr is correct
#### If you need the rows for first event analysis, take evt1=1
length(unique(SNs_py$sjlid[SNs_py$event==1]))
length(unique(SNs_py$sjlid[SNs_py$event==1 & SNs_py$evt1==1 ]))

SNs_py$PY <- SNs_py$end-SNs_py$start

###############
## Model fit ##
###############
# cc <- final[final$evt1==1,]

PHENO.ANY_SN <- SNs_py

## Age attained age (cubic spline)
source("/attributable_fraction_R_codes/cubic_spline.r")
breaks = seq(5, 95, 22.5)
cp = quantile(PHENO.ANY_SN$end, breaks/100, na.rm = T)
cs = cubic_spline(PHENO.ANY_SN$end, knots = cp)
PHENO.ANY_SN <- PHENO.ANY_SN[!colnames(PHENO.ANY_SN) %in% c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")]
colnames(cs) <- c("AGE_AT_LAST_CONTACT.cs1", "AGE_AT_LAST_CONTACT.cs2", "AGE_AT_LAST_CONTACT.cs3", "AGE_AT_LAST_CONTACT.cs4")
PHENO.ANY_SN <- cbind.data.frame(PHENO.ANY_SN, cs)

rm(list = setdiff(ls(), c("cc", "PHENO.ANY_SN")))

## **
save.image("/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.Any_SNs.V20b.Rdata")
