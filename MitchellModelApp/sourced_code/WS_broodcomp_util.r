# Code herein written by L. Mitchell (USFWS). This file contains the functions
# used in the model and used in file functions_mitchell_model.R (J. DuBois
# 10-Aug-2015)
# ******************************************************************************

## Utility functions and constants for use with the files 
## WS_broodcomp_data.r, WS_broodcomp_growthmodel.r, 
## WS_broodcomp_run.r, and possibly others.
##


##########################################################

library(pscl)       # zeroinfl()
library(MASS)       # glm.nb()
library(ALKr)       # contains the kimura_chikuni() fcn on which 
                    # kimura_chikuni_custom() is based

###########################################################
## Vectors for easier plotting:

## Plot colors for the length-only data catch years:
catchYrCol <- c("black","red","green3","blue","cyan","magenta","yellow","gray",
                "brown1","darkorange2","darkturquoise","deeppink1",
                "hotpink","khaki","lightgreen","lightslateblue","midnightblue",
                "mediumorchid","royalblue2",
                "peachpuff","skyblue2","thistle3","springgreen","yellow2",
                "steelblue4")
names(catchYrCol) <- as.character(1990:2014)


## Plot colors for brood years. Too many brood years, so color according
## to decade.
broodYrCol.decade <- c("black","red","green3","blue","cyan","magenta","yellow")
names(broodYrCol.decade) <- c("5","6","7","8","9","0","1")


## Plot colors for age-length data based on source.
sourceCol <- c("CDFW"="black", "DERBY"="orange", "OtterTrawlKey"="blue")



###########################################################
## Vectors and functions related to growth model-fitting:

growthCurveFcn <- function(paramObject, age) {
  # param can be an nls object or a vector.
  
  if(class(paramObject) == "nls") {
    paramVec <- coef(fit)
  } else {
    paramVec <- paramObject
  }

  if("K" %in% names(coef)) {  # von Bertalanffy
    Linf <- paramVec["Linf"]
    K <- paramVec["K"]
    t0 <- paramVec["t0"]
    
    return( Linf*(1-exp(-K*(age-t0))) )
    
  } else if("alpha" %in% names(coef)) { # Gompertz
    Linf <- paramVec["Linf"]
    alpha <- paramVec["alpha"]
    beta <- paramVec["beta"]

    return( Linf*exp(alpha*exp(beta*age)) )
  }
}




fit.GM <- function(growthModelType, plotBool=FALSE, newDev=TRUE) {
  # For debugging: growthModelType <- "VB.Normal"
  
  ## Used globally: agelenData
   
  if(growthModelType == "VB.Normal") {
    ## VBGM, Additive errors:  
    gmFormula <- ForkLength ~ 280*(1-exp(-K*(Age-t0)))
    init <- list(K = 0.15, t0 = 0)  #list(Linf = 225, K = 0.15, t0 = 0) 
    
  } else if(growthModelType == "VB.LogNormal") {
    ## VBGM, Multiplicative errors:
    gmFormula <- ln_ForkLength ~ log( Linf*(1-exp(-K*(Age-t0))) )
    init <- list(Linf = 225, K = 0.15, t0 = -1)
       
  } else if(growthModelType == "Gompertz.Normal") {
    ## Gompertz, Additive errors:
    gmFormula <- ForkLength ~ Linf*exp(alpha*exp(beta*Age))
    init <- list(Linf = 225, alpha = -1.5, beta = -0.15) 
      
  } else if(growthModelType == "Gompertz.LogNormal") {
    ## Gompertz, Multiplicative errors:
    gmFormula <- ln_ForkLength ~ log( Linf*exp(alpha*exp(beta*Age)) )
    init <- list(Linf = 225, alpha = -1.5, beta = -0.15)
       
  }

  modelName <- strsplit(growthModelType, split=".", fixed=TRUE)[[1]][1]
  modelError <- strsplit(growthModelType, split=".", fixed=TRUE)[[1]][2]
       
  fit <- nls(gmFormula, data=agelenData, start=init)
  fit[["growthModelType"]] <- growthModelType
  fit[["modelName"]] <- modelName
  fit[["modelError"]] <- modelError
  
  ## Create plots:
  if(plotBool) {
    pred <- predict(fit)
    res <- resid(fit)
     
    if(newDev) dev.new(width=11, height=9, units="inch")
    par(mfrow=c(2,2),mar=c(4.5,4.5,2,2),oma=c(0,0,2,1),cex.axis=1.1,
          cex.lab=1.1)
    
    xLim <- c(0,70)
    yLim_N <- c(0,300)
    yLim_LN <- c(3,5.9) 
    if(modelError=="Normal") {
      plot(agelenData$Age, agelenData$ForkLength, xlab="Age (years)", 
            ylab="ForkLength (cm)", ylim=yLim_N, #max(agelenData$ForkLength)),
            col=sourceCol[agelenData$Source], xlim=xLim)
    } else if(modelError=="LogNormal") {    
      plot(agelenData$Age, agelenData$ln_ForkLength, xlab="Age (years)", 
            ylab="Log(ForkLength (cm))", ylim=yLim_LN, xlim=xLim,
            col=sourceCol[agelenData$Source])
    }
    legend("topleft", names(sourceCol), col=sourceCol, pch=1, title="Source")   

    lines(agelenData$Age, pred, col="red", lwd=2)
    lines(0:70, predict(fit,data.frame("Age"=0:70)), col="red", lwd=1, lty=2)
    hist(res, main="", xlab="Residuals")
    plot(pred, res, xlab="Predicted ForkLengths", ylab="Residuals")
    
    if(modelError=="LogNormal") {
      plot(agelenData$Age, agelenData$ForkLength, xlab="Age (years)", 
            ylab="ForkLength (cm)", ylim=yLim_N, #max(agelenData$ForkLength)),
            col=sourceCol[agelenData$Source], xlim=xLim,
            main="Original Scale")
      lines(agelenData$Age, exp(pred), col="green3", lwd=2)
      lines(0:70, exp(predict(fit,data.frame("Age"=0:70))), col="green3", lwd=1, lty=2)      
    }

    title(main="Age-Length Data", outer=TRUE, line=0.5)
    title(main=paste(modelName,"with",modelError,"errors"), outer=TRUE, line=-0.9)      
#    title(main="Age-Length Data", outer=TRUE, line=-2)
#    title(main=paste(modelName,"with",modelError,"errors"), outer=TRUE, line=-3.5)  
  }
  
  return(fit)
}


#plotCompareGM <- function() {
#  # Compare additive and multiplicative error:
#  dev.new(width=7, height=6, units="inch")
#  plot(or$finalage, or$forklengthcm, xlab="Age (years)", ylab="ForkLength (cm)"))
#
#  lines(or$finalage, pred_additive, col="blue", lwd=2)
#  lines(or$finalage, exp(pred_multi), col="green3", lwd=2)
#  
#  title(main="OR Age-Length Data", outer=TRUE, line=-2)
#  title(main="Comparison of VBGM curves", outer=TRUE, line=-3.5)
#  legend("bottomright", c("Additive","Multiplicative"), col=c("blue","green3"), lwd=2)
#}
#





###########################################################
###########################################################
## Function to generate p_{L|A} using length and age classes
## specified in interval_table:

InvALK.fcn <- function(growthModelFit, len, age, interval_table) {
  # For debugging:
  # growthModelFit <- myFit; len <- lengthVec[1]; age <- ageVec; interval_table <- alk_lookup

  ## age can be a vector, but len should be a scalar.

  modelName <- growthModelFit[["modelName"]]
  modelError <- growthModelFit[["modelError"]]
  
  MN <- predict(growthModelFit, data.frame("Age"=as.numeric(age)))
  SD <- summary(growthModelFit)[["sigma"]] 
  
  # there should only be one match!
  index <- with(interval_table, which(lower <= len & len < upper))
  if(length(index) > 1) {
    cat("Error in the function InvALK()\n")
    return(NULL)
  }

  low <- interval_table[index,"lower"] 
  upp <- interval_table[index,"upper"]

  smallestLength <- interval_table[1,"lower"]
  largestLength <- interval_table[nrow(interval_table),"upper"]
  
  if(modelError == "Normal") {
    scaleBy <- pnorm(largestLength, mean=MN, sd=SD) - 
                    pnorm(smallestLength, mean=MN, sd=SD)
    probs <- (pnorm(upp, mean=MN, sd=SD) - pnorm(low, mean=MN, sd=SD))/scaleBy
    
  } else if(modelError == "LogNormal") {
    MN_adj <- MN - SD^2/2
    
    scaleBy <- plnorm(largestLength, meanlog=MN_adj, sd=SD) - 
                    plnorm(smallestLength, meanlog=MN_adj, sd=SD) 
    probs <- (plnorm(upp, meanlog=MN_adj, sd=SD) - plnorm(low, mean=MN_adj, sd=SD))/scaleBy


#    dev.new()
#    plot(1:300, pnorm(log(1:300),MN,SD))    
#    lines(1:300, plnorm(1:300, MN, SD), lwd=2, col="green3")
#    tmp <- rlnorm(2000, MN_adj, SD)
#    abline(v=median(tmp), col="red", lwd=2)
#    abline(v=mean(tmp), col="blue", lwd=2)
  } else {
    return(NULL)
  } 
  
  return(probs)
}




###########################################################
## Function to calculate iterated ALK using the Kimura-Chikuni algorithm.
## Adapted from the package ALKr.

kimura_chikuni_custom <- function(x, fi2, threshold = 1e-04, maxiter = 20000, 
                              age_classes = colnames(x), 
                              length_classes = rownames(x), name = "", 
                              description = "") 
{
    invAlk <- x # calc_invALK(calc_ALK(x), fi1)     # wants l|a, dimensions L X A
    pj2 <- rep(1/ncol(x), ncol(x))
    criterion <- 1
    iter <- 0
    
    while (criterion > threshold & iter < maxiter) {
        pj2.old <- pj2
        alk <- t(t(invAlk) * pj2)/rowSums(t(t(invAlk) * pj2)) 

        ############################
        ## Can have problems with dividing by 0 if no growth model is 
        ## used and there are holes in invALK.
        #cat("filling in nan's\n")
        alk[is.nan(alk)] <- 0
        ############################
                     
        nij <- fi2 * alk
        pj2 <- colSums(nij)/sum(nij)
        criterion <- sum(abs(pj2 - pj2.old))
        iter <- iter + 1
    }
    new("ALKr", alk = calc_ALK(nij), N = nij, method = "Kimura & Chikuni", 
        parameters = list(ConvergenceThreshold = threshold, Iterations = iter, 
            Converged = iter < maxiter), name = name, description = description)            
}





###########################################################


ageAssignResults <- function(lenDataSource, lengthBinVec, ageBinVec, alkType,
                             growthModelType="None") {
  ## Options for inputs:
  ##   lenDataSource: "Tagging", "Card", or "Both"
  ##   alkType: "Raw" or "Iterated"
  ##   growthModelType: "VB.Normal", "VB.LogNormal", "Gompertz.Normal", 
  ##                    "Gompertz.LogNormal", or "None"
  ## 
  ## lengthBinVec and ageBinVec have the structure:
  ##              c(min, practicalMax, absoluteMax, stepSize)
  ## 
  ## Examples:
  ##   lengthBinVec = c(20,250,300,5)
  ##   ageBinVec = c(0,40,70,1)
  ##  
  
  # For debugging:
#  lenDataSource <- "Tagging"; lengthBinVec=c(20,160,300,5); ageBinVec <- c(0,40,70,1)

#  alkType <- "Iterated"; growthModelType <- "None"
    
  # ****************************************************************************
  # changed below so lenDataSource input would accommodate more than the 3
  # selections # Validate (some) input:
  #if(!lenDataSource %in% c("Tagging","Card","Both")) {
  #  cat("lenDataSource must be 'Tagging', 'Card', or 'Both'\n")
  #  return(NULL)
  #} else if(!alkType %in% c("Raw","Iterated")) {
  #  cat("alkType must be 'Raw' or 'Iterated'\n")
  #  return(NULL)
  #} 
  # ****************************************************************************
  
  # ensure argument is supplied to lenDataSource
  if (missing(lenDataSource)) {
    stop("Argument missing. Supply dataframe or name of data source ", 
         "to 'lenDataSource'.", call. = FALSE)
  }
  
  # ensure argument passed to 'alkType' is either 'Raw' or 'Iterated'
  if (!(alkType %in% c("Raw", "Iterated"))) {
    stop("alkType argument must be 'Raw' or 'Iterated'", call. = FALSE)
  }
  
  validGM <- c("VB.Normal","VB.LogNormal","Gompertz.Normal","Gompertz.LogNormal","None")
  if(alkType=="Iterated" & !growthModelType %in% validGM)  {
    cat("Invalid growthModelType. See the function ageAssignResults().\n")
    return(NULL)              
  }
  
  ## Use the data frames agelenData cardData, taggingData, and bothData globally. 
  ## Define lenUse based on lenDataSource.
  agelenUse <- agelenData
  
  #if(lenDataSource == "Tagging") {
  #  lenUse <- taggingData
  #} else if(lenDataSource == "Card") {
  #  lenUse <- cardData
  #} else if(lenDataSource == "Both") {
  #  lenUse <- bothData
  #}
  
  # ******* added switch statement below in lieu of if/if-else above ***********
  # the switch is cleaner and offers option to use data from application in
  # GearSelectivityModelingApp.Rmd
  
  # needed to deparse(substitute()) lenDataSource if lenDataSource = data.frame
  # 'cause EXPR in switch does not accept objects length > 1
  data_source <- if (is.character(lenDataSource)) {
    lenDataSource
  } else {
    deparse(substitute(lenDataSource))
  }
  
  switch(EXPR = data_source,
         "Tagging"   = lenUse <- taggingData,
         "Card"      = lenUse <- cardData,
         "Both"      = lenUse <- bothData,
         "TagApp"    = lenUse <- taggingData(),
         "CardApp"   = lenUse <- cardData(),
         "BothApp"   = lenUse <- bothData(),
         lenUse <- lenDataSource # default - allows dataframe to be supplied
         )
  # ****************************************************************************

  
  ## Define length and age classes:
  lengthMin <- lengthBinVec[1]
  lengthMax <- lengthBinVec[2]    # left end of the last bin
  lengthPlus <- lengthBinVec[3]   # right end of the last bin; defines the "plus" length group
  lengthBinWidth <- lengthBinVec[4]
      
  ageMin <- ageBinVec[1]
  ageMax <- ageBinVec[2]          # left end of the last bin
  agePlus <- ageBinVec[3]         # right end of the last bin; defines the "plus" age group
  ageBinWidth <- ageBinVec[4]
    
  lengthBreaks <- c(seq(lengthMin, lengthMax, lengthBinWidth), lengthPlus)
  nLength <- length(lengthBreaks) - 1
  lengthVec <- (lengthBreaks[1:nLength] + tail(lengthBreaks,-1))/2
  lengthClass <- cut(lengthVec, breaks=lengthBreaks, right=FALSE)

  ageBreaks <- c(seq(ageMin, ageMax, ageBinWidth), agePlus)
  nAge <- length(ageBreaks) - 1
  ageVec <- ageBreaks[1:nAge]  # (ageBreaks[1:nAge] + tail(ageBreaks,-1))/2 
  ageClass <- cut(ageVec, breaks=ageBreaks, right=FALSE)
  

  ## Check if any observed ages or lengths are being excluded:  
  excludingLengths_1 <- with(lenUse, any(ForkLength < lengthMin) || 
                                            any(ForkLength > lengthPlus) )
  excludingLengths_2 <- with(agelenUse, any(ForkLength < lengthMin) || 
                                            any(ForkLength > lengthPlus) )                                            
  excludingAges <- with(agelenUse, any(Age < ageMin) || any(Age > agePlus) )
  
  if(excludingLengths_1) {
    warning("Choice of lengthBinVec is excluding length-only fish")
  }
                                          
  if(excludingLengths_2) {
    warning("Choice of lengthBinVec is excluding age-length fish")
  }
  
  if(excludingAges) {
    warning("Choice of ageBinVec is excluding age-length fish")
  }
  
    
  
  ## Define the length-only catch years and brood years:
  # Catch years
  year.catch <- sort(unique(lenUse$Year))
  year.catch.char <- as.character(year.catch)
  nYearCatch <- length(year.catch)
  
  # Brood years:
  year.brood <- sort(unique(c(sapply(year.catch, function(x) x - ageVec))))
  year.brood.char <- as.character(year.brood)
  nYearBrood <- length(year.brood)


  
  ## Bin agelenUse and lenUse according to length and age classes:
  agelenUse$ForkLength <- floor(agelenUse$ForkLength)
  agelenUse$lengthClass <- cut(agelenUse$ForkLength, breaks=lengthBreaks, right=FALSE)
  agelenUse$ageClass <- cut(agelenUse$Age, breaks=ageBreaks, right=FALSE)

  agelenBinned <- with(agelenUse, table(lengthClass, ageClass))


  ## Create two versions of lenBinned, one each in 'long' and 'wide' format:
  lenUse$lengthClass <- cut(lenUse$ForkLength, breaks=lengthBreaks, right=FALSE)
  
  # Long format: data.frame with columns = lengthClass, Year, Count
  lenBinnedLong <- expand.grid("lengthClass"=lengthClass, "Year"=year.catch)
  row.names(lenBinnedLong) <- with(lenBinnedLong, paste(Year, lengthClass, sep="."))
  lenBinnedLong[ ,"Count"] <- 0
  
  tmp <- with(lenUse, sapply(split(Count, list(Year,lengthClass), drop=TRUE), sum))  
  lenBinnedLong[names(tmp),"Count"] <- tmp
  
  # Wide format: matrix with rows = lengthClass and columns = catch year
  lenBinnedWide <- matrix(NA, nrow=nLength, ncol=nYearCatch, 
                             dimnames=list(lengthClass,year.catch.char))
  for(y.char in colnames(lenBinnedWide)) {
    tmp <- subset(lenBinnedLong, Year==as.numeric(y.char))
    tmpRowNames <- substring(row.names(tmp), first=6)
    lenBinnedWide[tmpRowNames,y.char] <- tmp$Count
  }

        
                             
  ## Construct ALK and apply it to get number-at-age and brood years.  
  nAgeList <- vector(mode="list", length=nYearCatch)
  names(nAgeList) <- year.catch.char
  broodList <- nAgeList
  
  if(alkType == "Raw") {
  
    # Use agelenBinned to directly construct an ALK, p_{A|L}:  
    alk <- matrix(0, nrow=nLength, ncol=nAge, dimnames=list(lengthClass,ageClass))
    for(L.char in rownames(agelenBinned)) {
      for(A.char in colnames(agelenBinned)) {
 
        # lengthClass and ageClass should be factors, because cut() was
        # used, in which case there won't be a problem here. 
        # Just to be sure, double check that A.char is represented in alk:     
        if(A.char %in% colnames(alk)) {
          alk[L.char, A.char] <- agelenBinned[L.char , A.char]
        }
        
      }
    }
    for(i in 1:nrow(alk)) {
      alk[i, ] <- alk[i, ] / sum(alk[i, ])  # Divide by the row total
    }
    alk[is.nan(alk)] <- 0                   # Replace NaNs with 0

    # Save the ALK:
    alkList <- list(alk)
        
    # Apply the ALK to each catch year:
    for(y in year.catch) {
      nLY <- lenBinnedWide[ ,as.character(y)]     
      nAY <- c(t(alk) %*% nLY)   # need [age x length]
      names(nAY) <- ageVec
      nAgeList[[as.character(y)]] <- nAY   

      tmp <- nAY
      names(tmp) <- (y - as.numeric(names(nAY)))
      broodList[[as.character(y)]] <- tmp  
    }
  
  } else if(alkType == "Iterated") {
  
    if(growthModelType != "None") {
      # Fit growthmodel:
      myFit <- fit.GM(growthModelType, plotBool=FALSE) 
      
      # This is used to integrate over length-at-age in invalk().
      alk_lookup <- data.frame("lower"=lengthBreaks[1:nLength],
                               "upper"=tail(lengthBreaks,-1))
  
      # Construct inverse ALK:
      invalk <- matrix(NA, nrow=nLength, ncol=nAge, dimnames=list(lengthClass,ageClass))
      for(i in seq(nrow(invalk))) {
        invalk[i, ] <- InvALK.fcn(myFit, lengthVec[i], ageVec, alk_lookup)
      }
      
    } else if(growthModelType == "None") {
      # Construct inverse ALK from the age-length data directly. There will be holes!
      invalk <- matrix(0, nrow=nLength, ncol=nAge, dimnames=list(lengthClass,ageClass))
      for(L.char in rownames(agelenBinned)) {
        for(A.char in colnames(agelenBinned)) {
        
          # Double check:
          if(A.char %in% colnames(invalk)) {
            invalk[L.char, A.char] <- agelenBinned[L.char , A.char]
          }
          
        }
      }
      for(i in 1:ncol(invalk)) {
        invalk[ ,i] <- invalk[ ,i] / sum(invalk[ ,i])   # Divide by the column total
      }
      invalk[is.nan(invalk)] <- 0                       # Replace NaNs with 0      
    
    }
  
    # Use inverse ALK to construct year-specific iterated ALKs: 
    alkList <- nAgeList
    
    for(y in year.catch) {
      nLY <- lenBinnedWide[ ,as.character(y)]
      ialk <- kimura_chikuni_custom(x=invalk, fi2=nLY)

      # Save the iterated ALK:
      alkList[[as.character(y)]] <- ialk
      
      nAY <- c(t(ialk@alk) %*% nLY)   # need [age x length]
      names(nAY) <- ageVec
      nAgeList[[as.character(y)]] <- nAY   
      # kc$N = pAgivenL_y with each column multiplied by nLY.
      # apply(kc@N,2,sum) is another way to get nAY.
      
      tmp <- nAY
      names(tmp) <- (y - as.numeric(names(nAY)))
      broodList[[as.character(y)]] <- tmp  
    }
                  
  } 
  

  ## Brood year counts by Brood Year and Catch Year (matrix):
  fullyrbrood <- unlist(lapply(broodList, function(x) as.numeric(names(x))))
  yrbrood <- sort(unique(fullyrbrood))   # should be identical to year.brood.char
  nBYxCY <- matrix(0, nrow=length(yrbrood), ncol=nYearCatch, 
                        dimnames=list(yrbrood,year.catch.char))
  for(y.char in colnames(nBYxCY)) {
    nBYxCY[names(broodList[[y.char]]),y.char] <- broodList[[y.char]]
  }
  
  ## Brood year counts by Brood Year only (vector): 
  nBY <- apply(nBYxCY, 1, sum)    
    

######################
######################
## This code will break if age classes aren't integers 0, 1, 2, ...

  ## Brood year counts by Brood Year and Age (matrix):
  cyear.byear <- unlist(broodList)
  split_by_by <- split(cyear.byear, substr(names(cyear.byear),6,9))
  split_by_by <- lapply(split_by_by, function(x) {
                    names(x) <- substr(names(x),1,4)
                    return(x) })
  for(i in names(split_by_by)) {
    tmp <- split_by_by[[i]]
    names(tmp) <- as.numeric(names(tmp)) - as.numeric(i)
    split_by_by[[i]] <- tmp
  }
  
  nBYxAgeWide <- matrix(NA, nrow=nYearBrood, ncol=length(ageVec), 
                        dimnames=list(year.brood.char,ageVec)) 
  for(by.char in rownames(nBYxAgeWide)) {
    ## The broodYear-age combinations that are impossible, given the 
    ## range of catch years, are left as NAs.
    ## The vector names(split_by_by[[by.char]]) ensures that only 
    ## physically possible ages get filled in.
      
    nBYxAgeWide[by.char, names(split_by_by[[by.char]])] <- split_by_by[[by.char]]
  }  

  ## A version of nBYxAgeWide that is long instead of wide:
  nBYxAgeLong <- expand.grid("Age"=as.numeric(colnames(nBYxAgeWide)),
                             "BroodYear"=as.numeric(rownames(nBYxAgeWide))) 
  for(i in seq(nrow(nBYxAgeLong))) {
    age.char <- as.character(nBYxAgeLong[i,"Age"])
    broodyr.char <- as.character(nBYxAgeLong[i,"BroodYear"])
    nBYxAgeLong[i,"Number"] <- nBYxAgeWide[broodyr.char,age.char]
  } 

######################
######################
                 
  ## Return lots of things.
  returnList <- list("lenDataSource"=data_source,#lenDataSource, 
                    "lengthBinVec"=lengthBinVec,
                    "ageBinVec"=ageBinVec, 
                    "alkType"=alkType,
                    "growthModelType"=growthModelType,
                    "lengthClass"=lengthClass, 
                    "ageClass"=ageClass,
                    "year.catch"=year.catch, 
                    "year.brood"=year.brood,
                    "agelenBinned"=agelenBinned, 
                    "lenBinnedLong"=lenBinnedLong,
                    "lenBinnedWide"=lenBinnedWide,
                    "alkList"=alkList, 
                    "nAgeList"=nAgeList, 
                    "broodList"=broodList,
                    "nBYxCY"=nBYxCY, 
                    "nBY"=nBY,
                    "nBYxAgeWide"=nBYxAgeWide, 
                    "nBYxAgeLong"=nBYxAgeLong)
                  
  if(alkType == "Iterated" & growthModelType != "None") {
    returnList[["gmFit"]] <- myFit
    returnList[["alk_lookup"]] <- alk_lookup
    returnList[["invalk"]] <- invalk
  }
  
  if(alkType == "Iterated" & growthModelType == "None") {
    returnList[["invalk"]] <- invalk
  }
  
  return(returnList)          
}





###########################################################


## Fit an exponential model to the number-at-age by brood-year:

fitExpModel <- function(resultObject, bySubsetVec, formulaOpt, plotBool=TRUE,
                          usedevnew=TRUE, display_vals = FALSE) {
  # For debugging: resultObject <- myResults; bySubsetVec <- c(1980:2013);
  
  #                resultObject <- myResults; bySubsetVec <- c(1980:2013);
  
  ## resultObject is a list returned by the ageAssignResults() function.
  ## bySubsetVec is a vector containing the brood years that should be
  ##   used in the model.
  ## formulaOpt is an integer indicating what type of model should be used,
  ##   see below.

  ## formulaOpt:
  ##  1 = NB; separate N0S0 for each brood year; one lambda.
  ##  2 = NB; one N0S0 for all brood years; one lambda.
  ##  3 = ZINB; separate N0S0 for each brood year; one lambda;
  ##            separate pi for each brood year. 
  ##  4 = ZINB; separate N0S0 for each brood year; one lambda;
  ##            one pi for all brood years.
  ##  5 = Poisson; separate N0S0 for each brood year; one lambda.
  ##  6 = ZIP; separate N0S0 for each brood year, one lambda; one pi.
  
  
  # Round the number-at-age to the nearest integer and create
  # a column with BroodYear as a factor:  
  nBYxAgeLong <- resultObject[["nBYxAgeLong"]]
  nBYxAgeLong$Num.int <- round(nBYxAgeLong$Number)
  nBYxAgeLong$BY.factor <- as.factor(nBYxAgeLong$BroodYear)
  
  # J. DuBois added (08-Jul-2015) because line 727:
  # barplot(round(t(nBYxAgeWide[by.char, ])) would not work without myResults
  # being run (from _run.r)
  nBYxAgeWide <- resultObject[["nBYxAgeWide"]]
  
  nBYxAgeLong_sub <- subset(nBYxAgeLong, BroodYear %in% bySubsetVec)
  
  
  # Define the model formula and fit it:
  if(formulaOpt == 1) {

    ## NB.
    ## Num.int_{BY} = N0S0_{BY} * exp(-Age*t)
    use.fmla <- as.formula(Num.int ~ BY.factor + Age)
    use.fit <- glm.nb(as.formula(use.fmla), data=nBYxAgeLong_sub,
                        na.action = na.omit) 
#                              control=glm.control(maxit=200,trace = 3))
                                
  } else if(formulaOpt == 2) {

    ## NB.
    ## Num.int_{BY} = N0S0 * exp(-Age*t)  
    use.fmla <- as.formula(Num.int ~ Age)
    use.fit <- glm.nb(as.formula(use.fmla), data=nBYxAgeLong_sub,
                        na.action = na.omit)
                                
  } else if(formulaOpt == 3) {
  
    ## Zero-inflated NB.
    ## Num.int_{BY} = N0S0_{BY} * exp(-Age*t)
    ## Separate pi for each brood year.    
    use.fmla <- as.formula(Num.int ~ BY.factor + Age | BY.factor)
    use.fit <- zeroinfl(use.fmla, dist="negbin", data=nBYxAgeLong_sub,
                        na.action = na.omit)
    
  } else if(formulaOpt == 4) {

    ## Zero-inflated NB.
    ## Num.int_{BY} = N0S0_{BY} * exp(-Age*t)
    ## One pi for all brood years.   
    use.fmla <- as.formula(Num.int ~ BY.factor + Age | 1)
    use.fit <- zeroinfl(use.fmla, dist="negbin", data=nBYxAgeLong_sub,
                        na.action = na.omit)
    
  } else if(formulaOpt == 5) {
  
    ## Poisson.
    ## Num.int_{BY} = N0S0_{BY} * exp(-Age*t)
    use.fmla <- as.formula(Num.int ~ BY.factor + Age)
    use.fit <- glm(as.formula(use.fmla), data=nBYxAgeLong_sub, family=poisson(),
                        na.action = na.omit) 
  
  } else if(formulaOpt == 6) {
  
    ## Zero-inflated Poisson.
    ## Num.int_{BY} = N0S0_{BY} * exp(-Age*t)
    ## One pi for all brood years.   
    use.fmla <- as.formula(Num.int ~ BY.factor + Age | 1)
    use.fit <- zeroinfl(use.fmla, dist="poisson", data=nBYxAgeLong_sub,
                        na.action = na.omit)  
      
  }


  # Separate N0S0 for each brood year:
  if(formulaOpt %in% c(1,3,4,5,6)) {
    # Unique, sorted being used to fit the model:
    uniq.by.char <- unique(as.character(nBYxAgeLong_sub$BroodYear)) 
         
    if(formulaOpt %in% c(1,5)) {   # non zero-inflated
      main_part <- use.fit$coef
    } else {                  # zero-inflated
      main_part <- use.fit$coef[["count"]]    # main parameters
      zero_part <- use.fit$coef[["zero"]]     # pi0-related parameters
    }

    lambda <-  -main_part[["Age"]]                    # the parameter called "Age"
    tmp <- main_part[!grepl("Age",names(main_part))]  # all other params
    N0S0 <- exp( c(tmp[1], tmp[1] + tmp[-1]) ) 
    names(N0S0) <- uniq.by.char

    if(plotBool) {
      plot_count <- 1
      
      # added par_old to reset default parameters in order to plot ggplots when
      # knitting to HTML (J.DuBois, 17-Jul-2015)
      par_old <- par(no.readonly = TRUE)
      
      for(by.char in uniq.by.char) {
        if(plot_count %% 9 == 1) {
          if(usedevnew) dev.new(width=11, height=7, units="inch")
          par(mfcol=c(3,3),mar=c(4.5,5,3,2),oma=c(0,0,2,1),cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
        }
  
        barplot(round(t(nBYxAgeWide[by.char, ])), xlab="Age", ylab="Number", 
                xlim=c(0,50), ylim=c(0,max(N0S0)), 
                main=paste0("Brood Year = ",by.char,", N0S0 = ",round(N0S0[by.char],1),
                            ", lambda = ",round(lambda,3)) )    
        lines(0:50, N0S0[by.char]*exp( -lambda*0:50), col="red", lwd=2)
  
        if(plot_count %% 9 == 1) {
          title(main=paste("Source:",resultObject[["lenDataSource"]]), outer=-1, 
                cex.main=1.5)
        }
           
        plot_count <- plot_count + 1
      }
      # added par(par_old) to reset default parameters in order to plot ggplots
      # when knitting to HTML (J.DuBois, 17-Jul-2015)
      par(par_old, pin = c(8, 8))
    }
  }

  # The same N0S0 for each brood year:
  if(formulaOpt == 2) {
    if(usedevnew) dev.new(width=11, height=6, units="inch")
    par(mfcol=c(1,1),mar=c(4.5,5,3,2),oma=c(0,0,2,1),cex.axis=1.3,cex.lab=1.3,cex.main=1.3)

    lambda <-  -use.fit$coef[["Age"]]
    N0S0 <- exp(use.fit$coef[["(Intercept)"]])
    
#    if(plotBool) {
#      plot(nBYxAgeLong_sub$Age, nBYxAgeLong_sub$Num.int,  
#              xlab="Age", ylab="Number", xlim=c(0,50), 
#              ylim=c(0,max(nBYxAgeLong_sub$Num.int)),
#              main=paste("Source:",resultObject[["lenDataSource"]]) )    
#      lines(0:50, N0S0*exp( -lambda*0:50), col="red", lwd=2)
#    }
  }


  # Barplot of just N0S0:
#  if(plotBool) {
#    if(usedevnew) dev.new(width=11, height=7, units="inch")
#    par(mfcol=c(1,1),mar=c(4.5,5,3,2),oma=c(0,0,2,1),cex.axis=1.3,cex.lab=1.3,cex.main=1.3)  
#    barplot(N0S0, xlab="Brood Year", ylab="Modeled Recruitment Index, N0S0")
#    title(main=paste("Source:",resultObject[["lenDataSource"]]), outer=-1, cex=1.5)
#  }    

  # Save values:
  use.fit[["aic"]] <- AIC(use.fit)         
  use.fit[["N0S0"]] <- N0S0
  use.fit[["lambda"]] <- lambda
      
  # J. DuBois added 28-Jul-2015 for more control over output
  if (display_vals) {
    # Display values:      
    cat("N0S0:\n");print(N0S0)
    cat("\nlambda:",lambda)
    cat("\nAIC:",use.fit[["aic"]],"\n")
  }
  
  return(use.fit)
}





wyPlots <- function(nBY_use, years=year.brood, 
                    xLab="Sacramento Water Year Index", 
                    yLab="Number caught by cohort",
                    mainOuter="") {
  # years can be numeric or character.
  years <- as.numeric(years)
  years.char <- as.character(years)
  x <- Sac_H20_Index[years.char]
  y <- nBY_use[years.char]
  Main <- paste(years[1],":",years[length(years)],sep="")
  scatter.smooth(x, y, pch=NA, xlab=xLab, ylab=yLab)
  text(x, y, paste0("'",substr(years.char,3,4)), cex=1.2 )
  title(main=Main, line=0.5)
  #title(main=mainOuter, outer=TRUE, line= -2)
}



myPlot <- function(x, y, ..., Height, Width, devNew=FALSE) {
  if(devNew) {
    dev.new(width=Width, height=Height, units="inch")
  }
  par(mfrow=c(3,3),mar=c(4.5,5,3,2),oma=c(0,0,2,1),cex.axis=1.3,cex.lab=1.3,
      cex.main=1.3)
  plot(x,y, ...)
}







#vb <- function(Linf, K, t0, Age) Linf*(1-exp(-K*(Age-t0)))
#
#gmFormula <- ForkLength ~ Linf*(1-exp(-K*(Age-t0)))
#init <- list(Linf=175.22, K=0.05, t0=-2.939) 
#tmpdata <- data.frame("Age"=0:50, "ForkLength"=
#            vb(init[[1]], init[[2]], init[[3]], 0:50) + rnorm(51,0,5))
#
#fit <- nls(gmFormula, data=tmpdata, start=init)
#fit[["growthModelType"]] <- "Custom.Normal"
#fit[["modelName"]] <- "Custom"
#fit[["modelError"]] <- "Normal"
#  





