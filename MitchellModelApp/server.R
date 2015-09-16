# This file contains the server for the MitchellModelApp (~July 2015)

# Libraries ---------------------------------------------------------------

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

# shinyServer -------------------------------------------------------------

# server is a work in progress as we make additional changes to the model

shinyServer(function(input, output) {
  
  # establish tagging dataframe using input from user-interface
  tagging_data <- reactive({
    GetStuTagGearSelData(
      catch_col = input$catch, millar_model = input$millar,
      min_len = input$minlen, l_mit_rows = TRUE)
  }) # end tagging_data()
  
  # establish card dataframe using input from user-interface
  card_data <- reactive({
    GetStuCardLFData(
      card_data = dfCardDataAll, legit_len = input$legitLen,
      species = "White", l_mit_rows = TRUE)
  }) # end card_data()
  
  agrument_list <- reactive({
    list(
      lenDataSource = input$data_source,
      lengthBinVec = c(0, 300, 310, 5),
      ageBinVec = c(0, 50, 70, 1),
      alkType = input$alk_type,
      growthModelType = input$growth_model
    )
  }) # end argument_list()
  
  results <- reactive({
    taggingData <- tagging_data()
    taggingData <<-
      subset(x = taggingData, subset = Year %in% input$years)
    cardData <<- card_data()
    
    GetAgeAssign(
      args_list = agrument_list(),
      fit_form_opt = input$fit_formula,
      show_plots = input$show_plt
    )
  }) # end results()
 
  #output$table <- renderTable({
    #testMe <- tagging_data()
    #tail(dfStuAll, 5)
    #agelenData # works
    #as.data.frame(Sac_H20_Index)# works too
    #head(testMe)
    #as.data.frame(agrument_list())
  #})
  
  output$summ <- renderPrint({
    cat("N0S0:\n")
    print(results()$N0S0)
    cat("\nLambda:\n", results()$lambda)
    cat("\nAIC:\n", results()$aic)
  })
  
  output$plot <- renderPlot({
    PlotN0S0(results()$N0S0)
  })
  
	
	
	
	plotByBroodYear_results <- reactive({
		plots <- plotByBroodYear(results(),3,3)		# nr and nc could be inputs; leave constant for now
		n <- length(plots)
		list(plots=plots, n=n)
	})

  output$plots_nByBroodYear <- renderUI({
    # do.call(tagList,
			lapply(1:plotByBroodYear_results()$n, function(i) {
				renderPlot({
					plotByBroodYear_results()$plots[[i]]
				})
				
		# )
		
    })
  })	
	
	
})
# end shinyServer
