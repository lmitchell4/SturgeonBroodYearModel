# Code herein written by L. Mitchell (USFWS) except where noted. This file
# contains the data for the 'artificial' BayStudy data and the age-length key
# data (J. DuBois 10-Aug-2015)
# ******************************************************************************

## Create data sets from the California White Sturgeon
## data files WST_Card.csv, WST_Tagging.csv, USFWSWSTDataAgeLen.csv.
## This file may also contain utility functions and constants.
##
## 
##

## agelenData   data frame containing California white sturgeon age-length data 
##                collected in 2014. The fish in the data set come from one of 
##                two sources, (1) CDFW or (2) the 2014 sturgeon derby (fish 
##                caught by participants and aged by Zac and team, I think). 
##                Each record represents one individual fish. The "Source" 
##                column indicates whether the fish came from CDFW of the DERBY.
##                Columns: Year, ForkLength, Age, BroodYear, ln_ForkLength, Source
## taggingData  data frame containing California white sturgeon length-only data
##                collected from routine trammel net surveys conducted by CDFW.
##                Each record represents a unique catchYear-forkLength category,
##                with the column "Count" indicating the number of fish in that 
##                category.
##                Columns: Year, ForkLength, Count
## cardData     data frame containing California white sturgeon length-only data
##                collected from angler report cards. Anglers are required to 
##                report the lengths of white sturgeon they keep, but not of 
##                white sturgeon they release, so report card lengths outside of
##                the legal slot are voluntarily reported.
##                Legal slot: 2007 to 2012(?): 46 - 66 inches total length,
##                            2013(?) to present: 40 - 66 inches fork length.
##                Has the same format as taggingData: each record represents a 
##                unique catchYear-forkLength category; the "Count" column 
##                indicates the number of fish in that category.
##                Columns: Year, ForkLength, Count
## bothData    data frame containing the California white sturgeon length-only
##                data from both the trammel net surveys and angler report cards. 
##                Has the same format as taggingData and cardData.
##                Columns: Year, ForkLength, Count


##########################################################


## Load age-length and length-only data:
#wsroot <- "C:/Main/Sturgeon/Sturgeon - White/Brood Year/"

# adding my directory in order to run code below using proper files (J. DuBois, 
# 30-Jun-2015); note: using same directory name for convenience
#wsroot <- "C:/Data/jdubois/RProjects/SturgeonBroodYearModel/DataFiles/"

# show files in directory
#dir(path = wsroot)

###### Artifical age-length points coming from the Bay Study Otter Trawl
###### age-month key. I think that key is based on total length, so 
###### convert to fork length using the formula Jason sent before.

#OtterAge <- c(rep(0,12), rep(1,12))

OtterAge <- c(1/12, 2/12, 3/12, 4/12, 5/12, 6/12, 7/12, 8/12, 9/12, 10/12, 11/12, 1,
           1 + 1/12, 1 + 2/12, 1 + 3/12, 1 + 4/12, 1 + 5/12, 1 + 6/12, 1 + 7/12,
           1 + 8/12, 1 + 9/12, 1 + 10/12, 1 + 11/12, 1 + 12/12)
OtterTotalLength <- c(8, 8, 8, 8, 16, 20, 24, 28, 32, 34, 36, 38,
                      38, 39, 40, 41, 42, 44, 46, 48, 50, 51, 52, 53)
OtterForkLength <- (0.9036 * OtterTotalLength) - 1.2162


BayStudyOtterTrawlKey <- data.frame(
  "Year" = NA,
  "ForkLength" = OtterForkLength,
  "Age" = OtterAge,
  "BroodYear" = NA,
  "ln_ForkLength" = log(OtterForkLength), # natural log
  "Source" = "OtterTrawlKey", 
  stringsAsFactors = FALSE
)




###### Clean 2014 Age-Length Data:
## Double check that the forklengths are in cm.
agelen_orig <- read.csv("data/USFWSWSTDataAgeLen.csv", header=TRUE,
                        stringsAsFactors=FALSE)

agelenData <- agelen_orig
agelenData$Year <- 2014      # as a reminder
agelenData$BroodYear <- with(agelenData, Year - Age)
agelenData$ln_ForkLength <- log(agelenData$ForkLength)

agelenData <- agelenData[ ,c("Year","ForkLength","Age","BroodYear",
                              "ln_ForkLength","Source")]

                              
## Combine agelenData and BayStudyOtterTrawlKey.
## If you don't want to use BayStudyOtterTrawlKey, just comment 
## the whole thing out above.
if(exists("BayStudyOtterTrawlKey")) {                
  agelenData <- rbind(BayStudyOtterTrawlKey, agelenData)
}

# Sort by Age, ForkLength, and Year:
agelenData <- agelenData[with(agelenData, order(Age,ForkLength,Year)), ]




###### Net Fathom Hour Effort Data:
nfhData <- read.csv("data/Sturgeon_AnnualEffortData.csv", header=TRUE,
										stringsAsFactors=FALSE)
row.names(nfhData) <- nfhData$Year







####### Water Year Index:
# Loads: Sac_H20_Index
#        SJ_H20_Index
#        Sac_H20_Index_adj
#        SJ_H20_Index_adj
#source(paste0("C:/Main/DSLCM/EDA/Water_Year_Data.R"))

# sourcing water year data from wsroot directory
#source(paste0(wsroot, "Water_Year_Data.R"))

# ******************************************************************************
# NOTE: in this copy of this file (WS_broodcomp_data2.r - a copy of 
# WS_broodcomp_data.r), I commented out the card and tagging data and function 
# containerFcn, as the card and tagging data will come from 
# GearSelectivityModelingApp.Rmd; did it this way 'cause this seemed the easiest
# at the time (J. DuBois 15-Jul-2015)
# ******************************************************************************

###### Clean Card Data:
## The original file contains all white sturgeon reported,
## regardless of whether or not the fork length was also reported.
## The legal slot limit is ~ 102 to 152 cm.
## Remove fish with fork length < 12 INCHES = 30.48 cm, per Jason's email.
## These are probably not accurate lengths (may accidentally be catches
## recorded in the wrong spot on the report card).
# card_orig <- read.csv(paste0(wsroot,"WST_Card.csv"), header=TRUE,
#                         stringsAsFactors=FALSE)
# 
# # Each record represents one fish:
# card_tmp <- subset(card_orig, !is.na(FL_cm) & FL_in >= 12)
# card_tmp$ForkLength <- card_tmp$FL_cm
# card_tmp$Count <- 1       # just a place holder for use in aggregate()
# 
# # Each record represents a different catchYear-forkLength:
# cardData <- aggregate(card_tmp[ ,"Count",drop=FALSE],
#                       by=card_tmp[ ,c("ForkLength","Year")],
#                       FUN=sum)
# cardData <- cardData[ ,c("Year","ForkLength","Count")]  # Reorder the columns
# row.names(cardData) <- with(cardData, paste(Year,ForkLength,sep="."))




###### Clean Tagging Data:
# tagging_orig <- read.csv(paste0(wsroot,"WST_Tagging.csv"), header=TRUE,
#                         stringsAsFactors=FALSE)
# 
# # Adding restraint on length because using gear selectivity model on small
# # lengths (where n is low) results in exceptionally high n after 'adjustment'
# # (J. DuBois, 06-Jul-2015) (working or not working? not sure J. DuBois)
# #tagging_orig <- tagging_orig[tagging_orig$TL_cm > 84, ]
# 
# tagging_orig$FL_cm <- (0.9036 * tagging_orig$TL_cm) - 1.2162    # from Jason's email
# 
# # tagging_orig has a separate count record for each unique catchYear-totalLength.
# # Round calculated forkLengths to the nearest integer, then re-tabulate based 
# # on unique catchYear-forkLength:
# tagging_tmp <- tagging_orig
# tagging_tmp$Count <- round(tagging_tmp$AdjCatch)
# tagging_tmp$ForkLength <- round(tagging_tmp$FL_cm)
# 
# # Each record represents a different catchYear-forkLength:
# taggingData <- aggregate(tagging_tmp[,"Count",drop=FALSE], 
#                           by=tagging_tmp[ ,c("ForkLength","Year")],
#                           FUN=sum)
# taggingData <- taggingData[ ,c("Year","ForkLength","Count")]  # Reorder the columns
# row.names(taggingData) <- with(taggingData, paste(Year,ForkLength,sep="."))




###### Create a data set containing both Card and Tagging Data:
# both_tmp <- rbind(cardData, taggingData)
# 
# # Each record represents a different catchYear-forkLength:
# bothData <- aggregate(both_tmp[,"Count",drop=FALSE], 
#                         by=both_tmp[ ,c("ForkLength","Year")],
#                         FUN=sum) 
# bothData <- bothData[ ,c("Year","ForkLength","Count")]  # Reorder the columns
# row.names(bothData) <- with(bothData, paste(Year,ForkLength,sep="."))








## Create a standardized verison of Sac_H20_Index:
#Sac_H20_Index.std <- (Sac_H20_Index - mean(Sac_H20_Index)) / sd(Sac_H20_Index)
#
#
#Sac_H20_Index_ma <- filter(Sac_H20_Index, c(1,2,0))    # moving average
#names(Sac_H20_Index_ma) <- names(Sac_H20_Index)
#
#
#
######## 
#source(paste0("C:/Main/DSLCM/Data/Clean_Data/PhysicalVariable_R_objects/",
#          "Entrain.Physical.Monthly.R"))
#
#
#meanOutflow.Mar.July <- sapply(split(Entrain.Physical.Monthly, 
#                                Entrain.Physical.Monthly$Year), 
#        function(x) {
#          tmp <- subset(x, Month %in% month.name[3:7])
#          return( mean(tmp$Outflow) )
#})
#



##########################################################






  
# containerFcn <- function() {
#   ## This is in here so that it can remain in this script, but
#   ## isn't run every time the script is sourced. I just don't
#   ## want to create a separate EDA script.
#   
#   ## EDA:
#   
#   plot_card_histogram <- function(savingToFile) {
#     ## Plot histogram of the Card fork length data.
#     
#     if(savingToFile) {
#       cexx <- 1.7
#       cexx.main <- 2.3
#     } else {
#       if(!savingToFile) dev.new(width=12, height=10, units="inch") 
#       cexx <- cexx.main <- 1.3
#     }
#   
#     par(mfrow=c(3,3),mar=c(4.5,4.5,3.5,2),oma=c(0,0,3,1),cex.axis=cexx,cex.lab=cexx,
#           cex.main=cexx.main)
#     xLim <- c(0, max(cardData$ForkLength))  # range(cardData$ForkLength)
#     yLim <- c(0,1100)
#     
#     for(y in unique(cardData$Year)) {
#       tmp <- subset(cardData, Year==y)
#       tmpVec <- with(tmp, rep(ForkLength, Count))
#       hist(tmpVec, xlab="ForkLength (cm)", main=y, xlim=xLim, ylim=yLim)
#       abline(v=c(102,152), col="red", lty=2)
#     } 
#     title(main="Card Data", outer=TRUE)
#     
#   }  
#     
#   
#   plot_tagging_histogram <- function(savingToFile) { 
#     ## Plot histogram of the Tagging fork length data.
#   
#     if(savingToFile) {
#       cexx <- 1.7
#       cexx.main <- 2.3
#     } else {
#       cexx <- cexx.main <- 1.3
#     }
#     
#     xLim <- c(0, max(taggingData$ForkLength) + 5)
#     #yLim <- c(0, max(taggingData$AdjCatch) + 1000)
#     yLim <- c(0,1000)
#     
#     count <- 1
#     for(y in unique(taggingData$Year)) {
#       if(count %% 9 == 1) {
#         if(!savingToFile) dev.new(width=12, height=10, units="inch")
#         par(mfrow=c(3,3),mar=c(4.5,4.5,3.5,2),oma=c(0,0,3,1),cex.axis=cexx,
#             cex.lab=cexx, cex.main=cexx.main)
#       }
#       tmp <- subset(taggingData, Year==y) 
#       tmpVec <- with(tmp, rep(ForkLength, Count))
#       hh <- hist(tmpVec, xlab="ForkLength (cm)", main=y, xlim=xLim, ylim=yLim)
#     
#       largeFreq <- hh[["counts"]][hh[["counts"]] > 1100]
#       if(length(largeFreq ) > 0) {
#         text(90, yLim[2]-50, paste0(largeFreq,collapse="     "), col="red", 
#               cex=cexx, adj=0)
#       }
#       
#       if(count %% 9 == 1) {
#         title(main="Tagging Data", outer=TRUE)
#       }
#           
#       count <- count + 1
#     } 
#   
#   }
# 
# 
# 
#   
#   ## Histogram of all length-only data:
#   dev.new(width=10, height=7, units="inch")
#   hist(with(bothData, rep(ForkLength,Count)), xlab="ForkLength (cm)",
#         main="Card and Tagging Data from All Years", xlim=c(0,350))
#         
#   agelenSplitBySource <- split(agelen_orig, agelen_orig$Source)
#   
#   dev.new(width=13, height=7, units="inch")
#   par(mfrow=c(1,2),mar=c(4.5,4.5,3,2),oma=c(0,0,1.5,1),cex.axis=1.2,cex.lab=1.2,
#         cex.main=1.2)
#   xLim <- range(agelen_orig$ForkLength)
#   yLim <- c(0,30)
#   plotFcn <- function(datTypeString, hatchDensity) {
#     with(agelenSplitBySource[[datTypeString]], hist(ForkLength, xlim=xLim, 
#           ylim=yLim, xlab="ForkLength (cm)", main=datTypeString))
#   }
#   plotFcn("CDFW")
#   plotFcn("DERBY")
#   title(main="2014", outer=TRUE)
#   
#   
#   
#   
#   
#   ###### Age-Length Data:
#   ## Plot length vs. age:
#   dev.new(width=10, height=7, units="inch")
#   with(agelenData, plot(Age, ForkLength, xlab="Age (years)", 
#         ylab="ForkLength (cm)", main="Age-Length Data"))
#         
#           
#   ###### Card Data:
#   plot_card_histogram(savingToFile=FALSE)
#   
#   ## Summary of Card Data (based on one record per fish):
#   summary(with(cardData, rep(ForkLength,Count)))
# 
# 
#   ###### Tagging Data:
#   plot_tagging_histogram(savingToFile=FALSE)
#   
#   ## Summary of Tagging Data (based on one record per fish):
#   summary(with(taggingData, rep(ForkLength,Count)))
# 
# 
# 
# 
# 
#   ###### Same plots as above, except save to pdf:
#   library(gplots)  
#   
#   pdf("WS.pdf", width=18, height=10, onefile=TRUE)
#   plot_card_histogram(savingToFile=TRUE)
#   
#   plot_tagging_histogram(savingToFile=TRUE)
#   
#   par(mfrow=c(1,1))
#     tmp1 <- as.matrix(t(summary(with(cardData, rep(ForkLength,Count)))))
#     rownames(tmp1) <- "Card"
#     
#     tmp2 <- as.matrix(t(summary(with(taggingData, rep(ForkLength,Count)))))
#     rownames(tmp2) <- "Tagging"
#     
#     tmp3 <- rbind(tmp1,tmp2)
#     textplot(tmp3, cex=2)  
#     
# #    textplot( as.matrix(t(summary(with(cardData, rep(ForkLength,Count))))) )
# #    textplot( as.matrix(t(summary(with(taggingData, rep(ForkLength,Count))))) )
# 
#   dev.off()
# 
# 
# 
# 
# 
# 
# }
# 
