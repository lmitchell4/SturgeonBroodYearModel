# 
# 9/16/15 (LM): Using this file for development purposes so I can run the
#   model in R and bypass Shiny. 
# ******************************************************************************


library(ggplot2)
library(plyr)
library(grid)
library(gridExtra)
library(gtable)
library(scales)
library(reshape)
library(shiny)
library(shinyapps)

# Source files and data ---------------------------------------------------

# data
load(file = "data/CDFWSturgeonData.RData")

# directories for source files 
source_dir <- "sourced_code/"
millar_files <- "millar_files/"

# source files
local_bool <- TRUE # sourcing files as local

# for application
source(file = paste0(source_dir, "functions_app_specific.R"),
       local = local_bool)

# for Mitchell Model
source(file = paste0(source_dir, "WS_broodcomp_util.r"),
       local = local_bool)
source(file = paste0(source_dir, "WS_broodcomp_data2.r"),
       local = local_bool)
source(file = paste0(source_dir, "Water_Year_Data.R"),
       local = local_bool)
source(file = paste0(source_dir, "functions_mitchell_model.R"),
       local = local_bool)

# for gear selectivity
source(file = paste0(source_dir, "functions_len_freq.R"),
       local = local_bool)
source(file = paste0(source_dir, "functions_gear_selectivity.R"),
       local = local_bool)

# for Millar functions used in gear selectivity
source(file = paste0(millar_files, "NextGeneration.R"),
       local = local_bool)
source(file = paste0(millar_files, "SelnCurveDefinitions.R"),
       local = local_bool)

# clean up
rm(source_dir, millar_files, local_bool)



# Select model input ---------------------------------------------------

# establish tagging dataframe using input from user-interface
taggingData <- GetStuTagGearSelData(
	catch_col = "AdjCatch",				# c("Catch", "AdjCatch")
	millar_model = "bilognorm",		# c("norm.loc", "norm.sca", "lognorm", "binorm.sca", "bilognorm")
	min_len = 40, 
	l_mit_rows = TRUE
)


# establish card dataframe using input from user-interface
card_data <- GetStuCardLFData(
		card_data = dfCardDataAll, 
		legit_len = 20,
		species = "White", 
		l_mit_rows = TRUE
)


args_list <- list(
	lenDataSource = "Tagging",					# c("Tagging", "Card")
	lengthBinVec = c(0, 300, 310, 5),
	ageBinVec = c(0, 50, 70, 1),
	alkType = "Iterated",				# c("Raw", "Iterated")
	growthModelType = "None"		# c("None", "VB.Normal", "VB.LogNormal", "Gompertz.Normal", "Gompertz.LogNormal")
)


results <- ageAssignResults(
	lenDataSource = args_list$lenDataSource,
	lengthBinVec = args_list$lengthBinVec,
	ageBinVec = args_list$ageBinVec,
	alkType = args_list$alkType,
	growthModelType = args_list$growthModelType
)


fit_exp <- fitExpModel(
	resultObject = results, 
	bySubsetVec = c(1970:2013),
	formulaOpt = 7, 
	usedevnew = FALSE,
	plotBool = FALSE
)


dev.new()
PlotN0S0(fit_exp$N0S0)




# results <- GetAgeAssign(
	# args_list = agrument_list,
	# fit_form_opt = 4,
	# show_plots = FALSE
# )

# PlotN0S0(results$N0S0)



