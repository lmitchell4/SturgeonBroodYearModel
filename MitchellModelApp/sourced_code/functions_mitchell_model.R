# Functions herein are used with the Lara Mitchell (LM, US Fish and Wildlife 
# Service) model for estimating a recruitment index of young White Sturgeon; I 
# wrote the functions herein in order to allow for more flexibility when 
# changing variables in the LM model; most of the functions are plotting
# functions (J. DuBois 08-Jul-2015)
# ******************************************************************************

GetAgeAssign <- function(args_list, fit_form_opt, wyi = Sac_H20_Index,
                         show_plots = FALSE) {
  # This function is a wrapper for the Lara Mitchell (LM) created
  # ageAssignResults(); this function basically does the same thing as the LM
  # creation but allows for easier passing of different arguments and plotting
  # (J. DuBois 07-Jul-2015)
  
  # Args:
  #   args_list: list of named arguments needed for the ageAssignResults()
  #              function (decided this would be the easiest way to pass 
  #              arguments to the ageAssignResults funtion)
  #
  #   fit_form_opt: 
  
  # Returns:
  #   not sure just yet (plots likely and some output of ageAssignResults(),
  #   which returns a list)
  
  # verify args_list is a list, stop otherwise
  if (!is.list(args_list)) {
    stop("'args_list' argument must be a list.", call. = FALSE)
  }
  
  # run function and assign output to results variable
  results <- ageAssignResults(
    lenDataSource = args_list$lenDataSource,
    lengthBinVec = args_list$lengthBinVec,
    ageBinVec = args_list$ageBinVec,
    alkType = args_list$alkType,
    growthModelType = args_list$growthModelType
  )
  
  # fit an exponential model to the number-at-age by brood year using the output
  # in the results variable
  fit_exp <- fitExpModel(
    # Lara used end year 2013 (should I automate this feature?); NOTE: using 
    # 2014 threw an error "Error in solve.default(as.matrix(fit$hessian)) system
    # is computationally singular: reciprocal condition number = 2.68247e-17" 
    # but I am not sure why (16-Jul-2015) ... has to do with zeroinfl(use.fmla, 
    # dist="negbin", data=nBYxAgeLong_sub, na.action = na.omit) in the 
    # "fitExpModel" function 'cause 2014 only has 1 row of data that is not NA
    # ... for now keep as 2013 will talk to L. Mitchell about this
    resultObject = results, bySubsetVec = c(1970:2013),
    formulaOpt = fit_form_opt, usedevnew = FALSE,
    plotBool = show_plots # FALSE
  )
  
	for(nm in names(fit_exp)) {
		results[[nm]] <- fit_exp[[nm]]
	}
	
  # display plots if desired *****************************************
  if (show_plots) {
    
    PlotBroodYr(by_x_cy = results$nBYxCY) # should yield two plots
    
    PlotN0S0(fit_exp$N0S0)
    
    PlotN0S0Wy(
      noso = fit_exp$N0S0, water_index = wyi,
      xLab = deparse(substitute(wyi)),
      yLab = "Modeled recruitment index, N0S0"
    )
    
    PlotAgeDist(results$nAgeList)
    
    # plot histogram of ages if age-length data exists in global environment
    if (exists("agelenData", envir = .GlobalEnv)) {
      
      # plot histogram of ages (plot is good for now but could set breaks in
      # geom_hist if desired)
      ggplot(data = agelenData, mapping = aes(x = Age)) +
        geom_histogram()
      
    }
    
    PlotAgePerBY(results$nBYxAgeWide)
    
  }
  # end show_plots ***************************************************
  
  # produces a dataframe of catch year (CY), number caught (per CY), number aged
  # (per CY), and the difference (former - latter), plus displays scatter of 
  # latter as a function of former with lin reg line of slope = 1 and intercept 
  # = 0 (J. DuBois, 08-Jul-2015)
  cat("\n") # just adding space for neater output
  #GetLengthsAgedStats(data_tag = taggingData, by_x_cy = results$nBYxCY)
  
  # function output
  # TBD
  #results
  #fit_exp
		
	return( results )
	
}
# end GetAgeAssign

PlotBroodYr <- function(by_x_cy) {
  # this function accepts a matrix of brood year by catch year and plots the
  # results as a graph for each catch year
  
  # Args:
  #   by_x_cy: a matrix of counts by brood year for each catch year
  
  # Returns:
  #   a plot for each catch year as count as a function of brood year (BY)
  
  # stop if by_x_cy is not a matrix
  if (!is.matrix(by_x_cy)) {
    stop("by_x_cy must be a matrix.", call. = FALSE)
  }
  
  # create dataframe for count of brood year (BY) by catch year (CY) and reshape
  # for plotting with ggplot
  dfByxCY <- as.data.frame(by_x_cy)
  dfByxCY$BY <- row.names(dfByxCY)
  dfByxCY$BYD <- AssignDecade(row.names(dfByxCY)) # brood year decade
  dfByxCY_melt <- reshape::melt.data.frame(
    dfByxCY, id.vars = c("BY", "BYD"),
    variable_name = "CY"
  )
  
  dfByxCY_melt$CY <- as.character(dfByxCY_melt$CY)
  
  # create range for plotting (x-axis breaks and minor breaks)
  by_range <- range(as.numeric(dfByxCY_melt$BY), na.rm = TRUE) # brood
  cy_range <- range(as.numeric(dfByxCY_melt$CY), na.rm = TRUE) # catch
  
  # create plot for each catch year of count as a function of brood year
  by_plot <- ggplot(data = dfByxCY_melt,
                    mapping = aes(x = as.numeric(BY), y = value)) +
    geom_bar(stat = "identity") +
    facet_wrap(facets = ~CY, nrow = 6) +
    scale_x_continuous(breaks =
                         seq(from = by_range[1], to = by_range[2], by = 10),
                       minor_breaks =
                         seq(from = by_range[1], to = by_range[2], by = 1),
                       expand = c(0, 0.5)) +
    labs(x = "Brood year", y = "Estimated count")
  
  # create barplot for count as a function of catch year (where bars are colored
  # according to decade)
  cy_plot <- ggplot(data = dfByxCY_melt,
                    mapping = aes(x = as.numeric(CY), y = value)) +
    geom_bar(mapping = aes(fill = BYD), stat = "identity") +
    #facet_wrap(facets = ~CY, nrow = 6) +
    scale_x_continuous(breaks =
                         seq(from = cy_range[1], to = cy_range[2], by = 2),
                       minor_breaks =
                         seq(from = cy_range[1], to = cy_range[2], by = 1),
                       expand = c(0, 0.5)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    labs(x = "Catch year", y = "Estimated count", fill = "BY Decade")
  
  # function output
  print(by_plot)
  print(cy_plot)
  #str(dfByxCY_melt)
}
# end PlotBroodYr

PlotN0S0 <- function(fit_exp_ns) {
  # This funtion plots (barplot) the modeled recruitment index (N0S0) by a
  # subsetted brood year
  
  # Args:
  #   fit_exp_ns: N0S0 from the fitExpModel() output
  
  # Returns:
  #   plot of modeled recruitment index as a function of brood year
  
  dfN0S0 <- data.frame(year = as.numeric(names(fit_exp_ns)),
                       mod_recruit = fit_exp_ns,
                       row.names = 1:length(fit_exp_ns),
                       stringsAsFactors = FALSE)
  
  # create range for plotting (x-axis breaks and minor breaks)
  by_range <- range(dfN0S0$year, na.rm = TRUE)
  
  # create plot to be displayed (bar plot)
  plot_ns <- ggplot(data = dfN0S0, mapping = aes(x = year, y = mod_recruit)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(breaks =
                         seq(from = by_range[1], to = by_range[2], by = 5),
                       minor_breaks =
                         seq(from = by_range[1], to = by_range[2], by = 1),
                       expand = c(0, 0.5)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    labs(x = "Brood year", y = "Modeled recruitment index, N0S0")
  
  # function output
  print(plot_ns)
}
# end PlotN0S0

PlotN0S0Wy <- function(noso, water_index, xLab = NULL, yLab = NULL) {
  # This function plots a scatter of N0S0 as a function of water year index
  
  # Args:
  #   noso: modeled recruitment index (N0S0) from exponential fit output
  #   water_index: the desired water year index data
  
  # Returns:
  #   a scatterplot with loess line (and 95% CI) of modeled recruitment index
  #   as a function of water year index (sac or san joaquin); 2-digit number
  #   displayed on datapoint is brood year (see bySubsetVec)
  
  # get labels for x and y axes
  x_lab <- if (is.null(xLab)) {
    toupper(deparse(substitute(water_index)))
  } else {
    xLab
  }
  
  y_lab <- if (is.null(yLab)) {
    deparse(substitute(noso))
  } else {
    yLab
  }
  
  # create dataframe for ease of plotting
  dfWYN0S0 <- data.frame(year = as.numeric(names(noso)),
                         noso = noso,
                         wyin = water_index[names(noso)],
                         row.names = 1:length(noso),
                         stringsAsFactors = FALSE
  )
  
  # create plot N0S0 as a function of water year index
  noso_wy_plot <- ggplot(data = dfWYN0S0, mapping = aes(x = wyin, y = noso)) +
    stat_smooth(method = "loess") +
    geom_point(shape = 21, fill = "white", colour = "white", size = 10) +
    geom_text(mapping = aes(label = substr(x = year, 3, 4)),
              vjust = 0.4, hjust = 0.5, size = 5, fontface = 2) +
    labs(x = x_lab, y = y_lab)
  
  # function output
  print(noso_wy_plot)
}
# end PlotN0S0Wy

PlotAgeDist <- function(age_list) {
  # This function plots estimated count as a function of age for each catch year
  
  # Args:
  #   age_list: list where each element is catch year and data is estimated
  #             count for each age (~0 to ~50, ages may vary with input)
  
  # Returns:
  #   plot of estimated count as a function of age for each catch year
  
  # ensure input is a list (could really accept a dataframe too and then change
  # ldply code accordingly, but for now this will do fine)
  if (!is.list(age_list)) {
    stop("Argument must be a list.", call. = FALSE)
  }
  
  # create dataframe from list
  dfAgeDist <- plyr::ldply(.data = age_list, .id = "CY")
  
  # reshape dataframe for ease of plotting
  dfAgeDist_melt <- reshape::melt.data.frame(
    dfAgeDist, id.vars = "CY",
    variable_name = "Age"
  )
  
  # age to integer for ease of plotting
  dfAgeDist_melt$Age <- as.integer(as.character(dfAgeDist_melt$Age))
  
  # create range for plotting (x-axis breaks and minor breaks)
  age_range <- range(dfAgeDist_melt$Age, na.rm = TRUE)
  
  # create plot as estimated count as a function of age
  age_plot <- ggplot(data = dfAgeDist_melt, mapping = aes(x = Age, y = value)) +
    geom_bar(stat = "identity", width = 0.75) +
    # need free_y for now but may remove later (J. DuBois 08-Jul-2015)
    facet_wrap(facets = ~CY, nrow = 6, scales = "free_y") +
    scale_x_continuous(breaks =
                         seq(from = age_range[1], to = age_range[2], by = 5),
                       minor_breaks =
                         seq(from = age_range[1], to = age_range[2], by = 1),
                       expand = c(0, 0.5)) +
    labs(x = "Age (years)", y = "Estimated count")
  
  # function output
  print(age_plot)
}
# end PlotAgeDist

PlotAgePerBY <- function(n_by_age) {
  # This function plots for each brood year where count > 0, count as a function
  # of age (J. DuBois 08-Jul-2015)
  
  # Args:
  #   n_by_age: a matrix of brood years (rows) by ages (columns) where data
  #             is count (number at each age per brood year)
  
  # Returns:
  #   for each brood year where count > 0, count as a function of age (removed
  #   plots where count == 0 so as not to clutter display with emplty plots);
  #   also returns vector of years for which count of age = 0 (nothing to plot)
  
  # use wide dataframe to see which rows have a count > 0
  has_count <- rowSums(n_by_age, na.rm = TRUE) != 0
  
  # create new dataset which only uses brood years where count > 0
  new_data <- as.data.frame(n_by_age[has_count, ])
  
  # make rownames (which are the brood years) a new field for ease of plotting
  new_data$BY <- rownames(new_data)
  
  # reshape new data for plotting
  new_data_melt <- reshape::melt.data.frame(
    new_data, id.vars = "BY",
    variable_name = "Age"
  )
  
  # change Age from factor to int for ease of plotting
  new_data_melt$Age <- as.integer(as.character(new_data_melt$Age))
  
  # create range for plotting (x-axis breaks and minor breaks)
  age_range <- range(new_data_melt$Age, na.rm = TRUE)
  
  # create plot
  age_by_plot <- ggplot(data = new_data_melt,
                        mapping = aes(x = Age, y = value)) +
    geom_bar(stat = "identity") +
    # fills by row (change BY to factor and then re-level to fill by column)
    facet_wrap(facets = ~BY, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks =
                         seq(from = age_range[1], to = age_range[2], by = 5),
                       minor_breaks =
                         seq(from = age_range[1], to = age_range[2], by = 1),
                       expand = c(0, 0.5)) +
    labs(x = "Age (years)", y = "Estimated count")
  
  # function output
  print(age_by_plot)
  
  # displays years where count of (number at) age is 0
  cat("\nBrood years (not plotted) as count of age = 0:\n")
  print(noquote(row.names(n_by_age[!has_count, ])))
}
# end PlotAgePerBY

GetLengthsAgedStats <- function(data_tag, by_x_cy) {
  # This function gets diagnostics about 'length-only' fish not aged because of
  # 'holes' in the age-length key; according to LM by using
  # BayStudyOtterTrawlKey the number not ages should be relatively low (J.
  # DuBois 08-Jul-2015)
  
  # Args:
  #   data_tag: tagging data (likely the dataframe taggingData)
  #   by_x_cy: matrix of counts brood year by catch year (likely nBYxCY)
  
  # Returns:
  #   list with number caught by catch year, number assigned an age by catch
  #   year, difference (former - latter), and boolean expressing if all
  #   catch years are covered (TRUE if yes)
  
  # get number caught per catch year
  num_caught_cy <- with(data = data_tag, expr = {
    tapply(X = Count, INDEX = list(Year), FUN = sum, na.rm = TRUE)
  })
  
  # get number assigned an age per catch year
  num_assigned_age_cy <- colSums(x = by_x_cy, na.rm = TRUE)
  
  # verify years are equivalent
  years_equal <- all(names(num_caught_cy) == names(num_assigned_age_cy))
  
  years <- if (years_equal) {
    names(num_caught_cy)
  } else {
    # not certain this will work as intended and has not been tested
    sort(unique(c(names(num_caught_cy), names(num_assigned_age_cy))))
  }
  
  # create dataframe from variables above
  df_lens <- data.frame(CatchYear = years, NumCaught = num_caught_cy,
                        NumAged = num_assigned_age_cy,
                        Diff = num_caught_cy - num_assigned_age_cy,
                        row.names = 1:length(years),
                        stringsAsFactors = FALSE)
  
  # create list of stats for output and understanding of number of fish (per LM)
  # 'excluded due to holes in the age-length data'
  out_list <- list(
    Data = df_lens,
    AllCatchYearsCovered = years_equal
  )
  
  # create plot of num_assigned_age_cy as a function of num_caught_cy
  plot_df <- ggplot(data = df_lens, mapping = aes(x = NumCaught, y = NumAged)) +
    geom_abline(intercept = 0, slope = 1, colour = "red", size = 2) +
    geom_point(shape = 21, fill = "white", colour = "white", size = 10) +
    geom_text(mapping = aes(label = substr(x = CatchYear, 3, 4)),
              vjust = 0.4, hjust = 0.5, size = 5, fontface = 2) +
    labs(x = "Number caught (CY)", y = "Number assigned age (CY)")
  
  # function output
  print(plot_df)
  out_list
}
# end GetLengthsAgedStats

GetModelStartingData <- function(catch_col = c("Catch", "AdjCatch"),
                                 min_len = 85, millar_model = "bilognorm",
                                 param_row_select = 3, l_mit_rows = TRUE,
                                 legit_len = 12, local_source = TRUE) {
  
  # This function was orginally called testFun() and a development copy can be
  # found in GearSelectivityModelingResults.R. I brought it here in order to
  # dynamically get the stugeon tagging data and sturgeon card data used in the
  # Mitchell Model; this version has slight changes over the devolpment copy (J.
  # DuBois 16-Jul-2015)
  
  # Args:
  #   catch_col:        gets either raw catch (Catch) or gear-selectivity
  #                     modeled catch (AdjCatch)
  #   min_len:          minimum length (total) in tagging data on which to
  #                     employ gear selectivity model
  #   millar_model:     the Millar gear selectivity model to use
  #   param_row_select: the row number for the parameters to use in the selected
  #                     gear selectivity model (make less confusing)
  #   l_mit_rows:       (TRUE) outputs dataframes of this function according to
  #                     L. Mitchell code with rownumbers as YYYY.len (year.length)
  #   legit_len:        Minimum (legitimate) length for angler-reported Card
  #                     data (sometimes anglers report 3 for length, which is
  #                     likely count not length in inches)
  #   local_source:     TRUE (default) sources all files to the environment
  #                     of this function
  
  # Returns:
  #   a list of dataframes (tagging, card, and both [a summary of both card and
  #   tagging])
  
  # file paths for needed source files
  source_files <- "C:/Data/jdubois/RSourcedCode/"
  project_files <- "C:/Data/jdubois/RProjects/GearSelectivity/"
  millar_files <- "MillarGearSelectivityCode/"
  
  # source *******************************
  source(file = paste0(source_files, "libraries_used_often.R"),
         local = local_source)
  source(file = paste0(source_files, "functions_data_import_export.R"),
         local = local_source)
  source(file = paste0(source_files, "SturgeonFunctions.R"),
         local = local_source)
  source(file = paste0(source_files, "GlobalRFunctions.r"),
         local = local_source)
  source(file = paste0(source_files, "SturgeonDataFrames.R"),
         local = local_source)
  source(file = paste0(source_files, "functions_len_freq.R"),
         local = local_source)
  source(file = paste0(source_files, "functions_gear_selectivity.R"),
         local = local_source)
  # file below formerly "SturgeonPopulationData.R"
  source(file = paste0(source_files, "data_connect_sturgeon_all.R"),
         local = local_source)
  
  # project ******************************
  source(file = paste0(project_files, millar_files, "NextGeneration.R"),
         local = local_source)
  source(file = paste0(project_files, millar_files, "SelnCurveDefinitions.R"),
         local = local_source)
  
  # clean up *****************************
  rm(source_files, project_files, millar_files)
  
  # for testing purposes *******************************************************
  #return(ls())
  #GetStuTagGearSelData(catch_col = "Adj", min_len = 0, l_mit_rows = T)
  #MeshCatch(stu_data = dfStuAll, year_end = 2014, min_len = 0,
  #          even_mesh_only = FALSE, plot_den = TRUE)
  # ****************************************************************************
  
  # get sturgeon tagging data
  data_tagg <- GetStuTagGearSelData(
    catch_col = catch_col, min_len = min_len, millar_model = millar_model,
    param_row_select = param_row_select, l_mit_rows = l_mit_rows,
    plot_n_data = FALSE 
    )
  
  # NOTE: in RelRetModelFit() plotting of Millar plots and display of Millar 
  # model output is supressed with plot_n_data = FALSE; change to TRUE if desire
  # to see Millar plots and Millar model fit parameters and deviance (J. DuBois,
  # 28-Jul-2015)
  
  # parameters and default arguments for GetStuTagGearSelData
  #(catch_col = c("Catch", "AdjCatch"),
  #min_len = 85, millar_model = "bilognorm",
  #param_row_select = 3, l_mit_rows = FALSE)
  
  # get White Sturgeon ("W") card data
  data_card <- GetStuCardLFData(
    card_data = dfCardDataAll, legit_len = legit_len, species = "W",
    l_mit_rows = l_mit_rows
    )
  #GetStuCardLFData(card_data = , legit_len = , species = , l_mit_rows = )
  
  # combine both tagging and card data (per Mitchell Model we need to sum catch
  # [Count] data after combining [rbind] the two dataframes)
  data_both <- plyr::ddply(
    .data = rbind(data_card, data_tagg), .variables = .(Year, ForkLength),
    .fun = summarise, Count = sum(Count, na.rm = TRUE)
    )
  
  # for testing **************
  #str(data_tagg)
  #str(data_card)
  #str(data_both)
  # **************************
  
  # function output (going this way for now, thought about assigment to global
  # env, but I think a list should work pretty well as these aren't big
  # dataframes - J. DuBois 16-Jul-2015)
  list(
    TaggingData = data_tagg,
    CardData    = data_card,
    TagCardData = data_both
  )
}
# end GetModelStartingData















plotByBroodYear <- function(modelResult, nr, nc) {
	# Creates plots of the estimated number caught by brood year and age. 
	# Each brood year is represented in a separate plot within a grid with nr rows
	# 	and nc columns. When the grid is full, a new device is opened.
	# 
  # Args:
	#   modelResult: a list returned by the function ageAssignResults().
  #   nr: Number of rows in the plot grid.
	#		nc: Number of columns in the plot grid.
  #
  # Returns:
  #  

	nBYxAgeLong <- modelResult[["nBYxAgeLong"]]
	selectSource <- modelResult[["lenDataSource"]]
	year.brood <- modelResult[["year.brood"]]
	nYearBrood <- length(year.brood)	
	

	ppp <- nr*nc	# plots per page

	numPage <- ceiling(nYearBrood/ppp)
	yearToPage <- rep(1:numPage, each=ppp)[1:nYearBrood]
	names(yearToPage) <- as.character(year.brood)
	nBYxAgeLong$plotPage <- yearToPage[as.character(nBYxAgeLong$BroodYear)]

	plots <- list()
	for(i in 1:numPage) {
		xmin <- min(nBYxAgeLong$Age)
		xmax <- max(nBYxAgeLong$Age)
		
		xsub <- subset(nBYxAgeLong, plotPage==i & !is.na(Number))	# to prevent ggplot warnings	
		xsub$BroodYear <- paste("Brood Year",xsub$BroodYear)
		
		## If there are fewer than ppp plots, add blank plots:
		n <- length(unique(xsub$BroodYear))
		if(n < ppp) {
			tmp <- sapply(1:(ppp - n), function(y) paste(rep(" ",y),collapse=""))
			useLevels <- c(sort(unique(xsub$BroodYear)), tmp)
			xsub$BroodYear <- factor(xsub$BroodYear, levels=useLevels, ordered=TRUE)
		}
		
		mygg <- ggplot(xsub, aes(x=Age, y=Number)) + geom_point() +
							xlim(xmin, xmax) + ylab("Number of Fish") + 
							theme(text = element_text(size = 15)) 
							#+ 
							#ggtitle(paste("Estimated Catch by Brood Year and Age\n","Source:",selectSource))	
							
		if(max(xsub$Number) == 0) {
			mygg <- mygg + ylim(0,1) + facet_wrap( ~ BroodYear, nrow=nr, ncol=nc, scales="free_y", drop=FALSE) 
		} else {
			mygg <- mygg + facet_wrap( ~ BroodYear, nrow=nr, ncol=nc, scales="free_y", drop=TRUE) 
		}

		plots[[i]] <- mygg
		#dev.new(width=11, height=6, units="inch")
		#print(mygg)
	}

	return(plots)
}															



