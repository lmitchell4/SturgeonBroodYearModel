# Code herein obtained from
# https://www.stat.auckland.ac.nz/~millar/selectware/RNext/

NetFit <- function(Data,Meshsize,x0,rtype="norm.loc",rel.power=NULL,
                   display_param = TRUE) {
	# Stop if Meshsize isn't sorted smallest to largest:
	if(sum(sort(Meshsize)==Meshsize)!=length(Meshsize))
		stop("Mesh size must be ascending order")
	if(is.null(rel.power)) rel.power=rep(1,length(Meshsize))
		Counts=as.data.frame(Data[,-1])  # I added as.data.frame()
	if(ncol(Counts)!=length(Meshsize))
		stop("Number of mesh sizes should be ",ncol(Counts))
		
	CountPropns <- Counts/apply(Counts,1,sum,na.rm=TRUE)
	fullfit.l <- sum(Counts*log(CountPropns),na.rm=TRUE)
	r <- selncurves(rtype) #Get selection curve function
	fit <- optim(x0,nllhood,Data=Data,Meshsize=Meshsize,r=r,rel.power=rel.power,
			hessian=T,control=list(trace=F))
	# 15-Jul-2015 (J. DuBois) added cleaner version of displaying paramters and
	# deviance and added option to show (with default = TRUE)
	if (display_param) {
	  cat("\n****************** Model fit parameters & deviance ***************\n")
	  cat("Model:", rtype, "|", "Parameters:", fit$par, "|",
	      "Deviance:", 2 * (fullfit.l + fit$value), "\n")
	  cat("******************************************************************\n")
	  #cat("Parameters=",fit$par,",    Deviance=",2*(fullfit.l+fit$value),"\n")
	}
	
	invisible(c(fit,deviance=deviance,rtype=rtype,rel.power=list(rel.power),
		   Meshsize=list(Meshsize),Data=list(Data))) 
}

nllhood <- function(theta,Data,Meshsize,r,rel.power) {
	lens <- Data[,1]
	Counts <- Data[,-1]
	rmatrix <- outer(lens,Meshsize,r,theta)
	rmatrix[is.na(Counts)] = NA #No fitted retention for missing meshsizes
	rmatrix <- t(t(rmatrix)*rel.power)
	phi <- rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
	nll <- -sum(Counts*log(phi),na.rm=TRUE)		# Using SELECT method (?)
	return(nll) 
}
  
Estimates <- function(fit) {
	require("msm")
	x <- fit$par
	varx <- solve(fit$hess) 
	names <- c("Mode(mesh1)","Std dev.(mesh1)")
	switch(fit$rtype,
		"norm.loc"={ pars=x; varpars=varx },
		"norm.sca"={ pars=x; varpars=varx },
		"lognorm"={
		  pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)))
		  varpars=deltamethod(list(~exp(x1-x2^2),
				   ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),x,varx,ses=F)},
		"binorm.sca"={
		  pars=c(x[1:4],exp(x[5])/(1+exp(x[5])))
		  names=c("Mode1(mesh1)","Std dev.1(mesh1)",
						"Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
		  varpars=deltamethod(list(~x1,~x2,~x3,~x4,~exp(x5)/(1+exp(x5))),
													x,varx,ses=F)},
		"bilognorm"={
		  pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)),
				 exp(x[3]-x[4]^2),sqrt(exp(2*x[3]+x[4]^2)*(exp(x[4]^2)-1)),
				 exp(x[5])/(1+exp(x[5])))
		  names=c("Mode1(mesh1)","Std dev.1(mesh1)",
						"Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
		  varpars=deltamethod(
			list(~exp(x1-x2^2),~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1)),
				 ~exp(x3-x4^2),~sqrt(exp(2*x3+x4^2)*(exp(x4^2)-1)),
				 ~exp(x5)/(1+exp(x5))),x,varx,ses=F)},  
		"tt.logistic"={
		  pars=c(-x[1]/x[2],2*(log(3))/x[2],exp(x[3])/(1+exp(x[3])))
		  names=c("L50","SR","p")
		  varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2,~exp(x3)/(1+exp(x3))),
								 x,varx,ses=F)},
		stop(paste("\n",fit$rtype, "not recognised, possible curve types are \n", 
			"\"norm.loc\", \"norm.sca\", \"lognorm\" \n", 
			"\"binorm.sca\", \"bilognorm\", and \"tt.logistic\"")) 
	)#End of switch
	
	estimates=cbind(pars,sqrt(diag(varpars)))
	colnames(estimates)=c("par","s.e.")
	rownames(estimates)=names
	return(estimates) 
}

PlotCurves <- function(fit,Meshsize=NULL,plotlens=NULL,standardize=TRUE,
                       show_plots = TRUE,...) {
	r <- selncurves(fit$rtype) #Get selection curve function
	if(is.null(plotlens)) plotlens=fit$Data[,1]
	if(is.null(Meshsize)) Meshsize=fit$Meshsize
	plot.title=switch(fit$rtype,
			"norm.loc"="Normal (common spread)",
			"norm.sca"="Normal",
			"lognorm"="Lognormal",
			"binorm.sca"="Bi-normal",
			"bilognorm"="Bi-lognormal",
			"tt.logistic"="Control and logistic","")
	rmatrix <- outer(plotlens,Meshsize,r,fit$par)
	rmatrix <- t(t(rmatrix)*fit$rel.power)
	if(standardize) rmatrix=rmatrix/max(rmatrix)
	# J. DuBois added 28-Jul-2015 to include option to supress plots if desired
	if (show_plots) {
	  # Added main  and legend in matplot (J. DuBois 19-Aug-2014)
	  matplot(plotlens,rmatrix,type="l",las=1,ylim=c(0,1),
	          xlab="Length (cm)",ylab="Relative retention", main=plot.title,...)
	  legend(x = 210, y = 0.75, col = 1:length(Meshsize), lty = 1:length(Meshsize), legend = Meshsize/2.54, xpd = TRUE,
	         cex = 0.75, horiz = FALSE, pt.cex = 0.75, xjust = -0.75, yjust = 0.25, ncol = 1)
	}
	
	#abline(h=seq(0,1,0.25),lty=3)
	lenrmatrix=cbind(plotlens,rmatrix)
	colnames(lenrmatrix)=c("Length",Meshsize)
	invisible(lenrmatrix) 
  #head(lenrmatrix)
}

Summary <- function(fit,label="Deviance residuals", show_plots = TRUE,
                 xlabel="Length (cm)",ylabel="Mesh size (cm)",cex=1) {
	r <- selncurves(fit$rtype) #Get selection curve function
	lens <- fit$Data[,1]; nlens=length(lens)
	Meshsize <- fit$Meshsize; nmeshes=length(Meshsize)
	O <- fit$Data[,-1]; #Matrix of observed counts
	rmatrix <- outer(lens,Meshsize,r,fit$par)
	rmatrix[is.na(O)] = NA #No fitted retention for missing meshsizes
	rmatrix <- t(t(rmatrix)*fit$rel.power)
	phi <- rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
	E <- apply(O,1,sum,na.rm=TRUE)*phi #Matrix of expected counts
	Pearson.resids <- (O-E)/sqrt(E)
	Pearson.chisq <- sum(Pearson.resids^2,na.rm=TRUE)
	
	
	# https://onlinecourses.science.psu.edu/stat504/node/82
	# http://pic.dhe.ibm.com/infocenter/spssstat/v20r0m0/index.jsp?topic=%2Fcom.ibm.spss.statistics.help%2Falg_genlog-poisson_residuals_deviance.htm
	# http://www.rni.helsinki.fi/~boh/Teaching/GLMs/glmsl11.pdf
	
	wk <- O*log(O/E); wk[is.na(wk)]=0
	Dev.resids <- sign(O-E)*sqrt(2*(E-O+wk))
	Deviance <- sum(Dev.resids^2,na.rm=TRUE)
	full.l <- sum(-O+O*log(O),na.rm=TRUE)
	null.E <- matrix(apply(O,1,mean,na.rm=TRUE),nrow=nlens,ncol=nmeshes)
	null.l <- sum(-null.E+O*log(null.E),na.rm=TRUE)
	model.l <- sum(-E+O*log(E),na.rm=TRUE)
	NonZeroDat <- O[apply(O,1,sum,na.rm=TRUE)>0,]
	d.o.f. <- nrow(NonZeroDat)*(nmeshes-1)-length(fit$par)-sum(is.na(NonZeroDat))
	out <- rbind(null.l,model.l,full.l,Deviance,Pearson.chisq,d.o.f.)
	AreLensUnique <- (length(lens)==length(unique(lens)))
	# J. DuBois added 28-Jul-2015 to include option to supress plots if desired
	if (show_plots) {
	  if(nmeshes>2&AreLensUnique) {
	    plot(1,1,xlim=range(lens),xlab=xlabel,ylab=ylabel,
	         ylim=range(Meshsize)+(cex/50)*c(-1,1)*(max(Meshsize)-min(Meshsize)),
	         yaxt="n",type="n",main=label)
	    axis(2,Meshsize,Meshsize,las=1)   
	    for(i in 1:nlens)
	      for(j in 1:nmeshes)
	        points(lens[i],Meshsize[j],pch=ifelse(Dev.resids[i,j]>0,16,1),
	               cex=3*abs(Dev.resids[i,j])*cex/(abs(max(Dev.resids)))) 
	  } else {
	    if(nmeshes==2) {
	      Dev.resids.len=sign(Dev.resids[,2])*sqrt(apply(Dev.resids^2,1,sum))
	      plot(lens,Dev.resids.len,type=ifelse(AreLensUnique,"h","p"),las=1,
	           main=label,xlab=xlabel,ylab=ylabel,cex=cex)
	      abline(h=0) 
	    }
	  }
	}
	
	return(out)
}

# Created function to compare effects of effort (equal or not equal) and to plot
# adjusted catch (adjusted using the values of relative retention); J.DuBois 15-Sep-2014
AdjCatchCompare <- function(catch.data, rel.ret.equal = Plot.equal, rel.ret.diff = Plot.diff, yr = NULL){
  # This function adjusts catch using relative retention values obtained using
  # sturgeon catch data from 1990:2013; catch.data is adjusted by dividing catch
  # by relative retention values
  
  # Gets range of test data just FYI
  lenRange <- range(catch.data$TL_c)
  
  # Gets data without lengths
  length.col <- catch.data[,"TL_c"]
  
  # Gets index values where length in test data is a match with length in 
  # rel.ret.equal data (assume rel.ret.equal lengths and rel.ret.diff lengths
  # the same); use all(Plot.diff[,1]==Plot.equal[,1]) to verify assumption
  lenMatch <- match(x = catch.data$TL_c, table = rel.ret.equal[,1])
  
  # Divides relative retention by catch to get adjusted catch
  adj.catch.equal <- catch.data[,-1] / rel.ret.equal[lenMatch,-1]
  adj.catch.diff <- catch.data[,-1] / rel.ret.diff[lenMatch,-1]
  
  # Gets proportions of adjusted catch (proportion by row (i.e., length))
  adj.catch.equal.prop <- prop.table(x = as.matrix(x = adj.catch.equal), margin = 1)
  adj.catch.diff.prop <- prop.table(x = as.matrix(x = adj.catch.diff), margin = 1)
  
  # Changes column names for merging
  colnames(adj.catch.equal.prop) <- paste(colnames(adj.catch.equal.prop), "e", sep = ".")
  colnames(adj.catch.diff.prop) <- paste(colnames(adj.catch.diff.prop), "d", sep = ".")
  
  # Creates data frames for analyses
  dfCatchAdj <- data.frame(Len = length.col, adj.catch.equal.prop, adj.catch.diff.prop)
  dfCatchAdjDiff <- data.frame(Len = length.col, adj.catch.equal - adj.catch.diff) # add .prop for proportions
  colnames(x = dfCatchAdjDiff) <- c("Len", "Mesh6", "Mesh7", "Mesh8")
   # creates data frame for plotting catch
  dfAllCatch <- data.frame(catch.data, adj.catch.equal, adj.catch.diff)
  colnames(dfAllCatch) <- c("Len",
                            # n=catch not adjusted, e=adusjted with equal effort, d=adjusted with unequal effort
                            paste0(rep(x = c("Mesh6_", "Mesh7_", "Mesh8_"), times = 3),
                                   rep(x = c("n","e","d"), each = 3)
                                   )
                            )
  
  # Melt data for plotting
  dfCatchAdj_melt <- melt(data = dfCatchAdj, id.vars = "Len")
  dfCatchAdjDiff_melt <- melt(data = dfCatchAdjDiff, id.vars = "Len")
  dfAllCatch_melt <- melt(data = dfAllCatch, id.vars = "Len")  
  
  # Create new variable names for plotting / facetting
  varname <- strsplit(x = as.character(x = dfCatchAdj_melt$variable), split = ".", fixed = TRUE)
  
  # Add new variables to dataframe
  dfCatchAdj_melt <- within(data = dfCatchAdj_melt, expr = {
    Mesh <- sapply(X = varname, FUN = "[", 1)
    AdjType <- sapply(X = varname, FUN = "[", 2)
  })
  
  dfAllCatch_melt <- within(data = dfAllCatch_melt, expr = {
    Mesh <- gsub(pattern = "\\D", replacement = "", x = variable)
    Effort <- gsub(pattern = "^Mesh..", replacement = "", x = variable)
  })
  
  # Create plots *********************************************************************************
  # plots adjust catch (proportion) comparing equal (e) and unequal (d) effort for each length by mesh
  pltBar <- ggplot(data = dfCatchAdj_melt, mapping = aes(x = Len, y = value)) +
    geom_bar(mapping = aes(fill = AdjType), stat = "identity", position = "dodge", width = 0.75) +
    facet_grid(facets = Mesh ~ .) +
    scale_fill_manual(name = "Effort", values = c(d = "#599ad3", e = "#f9a65a"),
                      labels = c(d = "Unequal", e = "Equal")) +
    xlab("Length (cm TL)") + ylab("Proportion of catch")
  
  # plots scatter plot of equal effort (e) vs unequal effort (d); adds line of slope 1 and int = 0
  pltScatter <- ggplot(data =  dfCatchAdj_melt) +
    geom_abline(intercept = 0, slope = 1, colour = "grey70", size = 2) +
    geom_point(mapping = aes(x = value[AdjType == "e"], y = value[AdjType == "d"],
                             colour = Mesh[AdjType == "e"]), alpha = 2/5) +
    guides(colour = guide_legend(override.aes = list(size = 4), title = NULL)) + #
    theme(legend.position = c(0.12,0.80)) +
    ylab("Unequal effort") + xlab("Equal effort")
  
  # plots the difference between adjusted catch (equal effort) and adjusted catch (unequal effort)
  pltDiff <- ggplot(data =  dfCatchAdjDiff_melt, mapping = aes(x = Len, y = value)) +
    geom_hline(yintercept = 0, colour = "grey", size = 1) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    #scale_y_continuous(limits = c(-0.10,0.10)) +
    facet_grid(facets = variable ~ .) +
    ylab("Proportion of catch (equal - unequal)") + xlab("Length (cm TL)")
  
  pltCatch <- ggplot(data = dfAllCatch_melt, mapping = aes(x = Len, y = value)) +
    geom_bar(stat = "identity", width = 0.75) + 
    #(mapping = aes(fill = Effort), shape = 21, alpha = 2/5) +
    facet_grid(facets = Effort ~ .) + #Mesh
    ylab("Catch") + xlab("Length (cm TL)") + ggtitle(yr)

  # Function output
#   cat("This plot represents the catch (as a proportion of total for each length) adjusted using
#       relative retention calculated from equal fishing effort and relative retention
#       calculated from unequal fishing effort.\n")
#   print(pltBar)
#   
#   cat("This scatter plot plots the catch (as a proportion of total for each length) adjusted using
#       relative retention calculated from equal fishing effort vs adjusted using relative retention
#       calculated from unequal fishing effort.\n")
#   print(pltScatter)
#   
#   cat("This plot shows the difference in catch (just catch, i.e., NOT as a proportion of total for each length)
#       adjusted using relative retention calculated from equal fishing effort and adjusted using relative retention
#       calculated from unequal fishing effort.\n")
#   print(pltDiff)
#   
#   cat("This plot shows the catch by mesh size adjusted using relative retention from unequal effort (d),
#       adjusted using relative retention from equal effort (e), and not adjusted (n).\n")
  print(pltCatch)
}

