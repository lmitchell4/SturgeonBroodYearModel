# This file contains the user-interface for the MitchellModelApp (~July 2015)


# shinyUI -----------------------------------------------------------------

# ui is a work in progress as we make additional changes to the model

# create variables for choices of some selectInputs
gm_choices <- c("None", "VB.Normal", "VB.LogNormal",
                "Gompertz.Normal", "Gompertz.LogNormal")

data_choices <- c("Tagging", "Card")#, "Both",
                  #"TagApp", "CardApp", "BothApp")

# create variable to hold names of Millar models for gear selectivity
millar_mod_choices <- c("norm.loc", "norm.sca", "lognorm",
                        "binorm.sca", "bilognorm")

tagging_years <- c(1990, 1991, 1993, 1994, 1997, 1998,
                   2001, 2002, 2005, 2006, 2007, 2008,
                   2009, 2010, 2011, 2012, 2013, 2014)

# create user-interface
shinyUI(fluidPage(
  #theme = "bootstrap.css",
  titlePanel("Mitchell Model App"),
  
  sidebarLayout(
    sidebarPanel(
      
      fluidRow(
        column(width = 12,
               # data source
               selectInput(
                 "data_source", label = "Data source:",
                 choices = data_choices,
                 selected = "Tagging", width = "125px"
               ),
               # age-length key option
               selectInput(
                 "alk_type", label = "Age-length key:",
                 choices = c("Raw", "Iterated"),
                 selected = "Iterated", width = "125px"
               ),
               # growth model option
               selectInput(
                 "growth_model", label = "Growth model:",
                 choices = gm_choices,
                 selected = "", width = "175px"
               )
        ) # end column
      ), # end fluidRow
      
      fluidRow(
        column(width = 12,
               # formula fit number option
               selectInput(
                 "fit_formula", label = "Model fit:",
                 choices = c(1:6),
                 selected = 4, width = "90px"
               ),
               # plotting option
               selectInput(
                 "show_plt", label = "Show plots:",
                 choices = c(TRUE, FALSE),
                 selected = FALSE, width = "90px"
               )
        ) # end column
      ),# end fluidRow
      
      fluidRow(
        column(width = 6,
               # catch
               selectInput(
                 "catch", label = "Catch data:",
                 choices = c("Catch", "AdjCatch"),
                 selected = "AdjCatch", width = "125px"
               ),
               # selectivity model
               selectInput(
                 "millar", label = "Millar model:",
                 choices = millar_mod_choices,
                 selected = "bilognorm", width = "125px"
               )
        ) # end column
      ), # end fluidRow
      
      fluidRow(
        column(width = 12,
               # sturgeon total length
               sliderInput(
                 "minlen", label = "Minimum length (cm TL):",
                 min = 40, max = 85, value = 40, step = 1
               ),
               # sturgeon total length
               sliderInput(
                 "legitLen", label = "Minimum length (inches TL):",
                 min = 10, max = 30, value = 20, step = 1
               )
        ) # end column
      ), # end fluidRow
      
      fluidRow(
        column(width = 12,
               checkboxGroupInput(
                 inputId = "years",
                 label = "Years:",
                 choices = tagging_years,
                 selected = tagging_years,
                 inline = TRUE
               )
        ) # end width
      ) # end fluidRow
      
    ), # end sidebarPanel
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("plot")), 
        tabPanel("Summary", verbatimTextOutput("summ")), 
        #tabPanel("Table", tableOutput("table")),
				tabPanel("More Plots", 
									h3("Estimated Catch by Brood Year and Age", align = "center"),
									uiOutput("plots_nByBroodYear")),
        selected = "Plot"
      ) # end tabsetPanel
      
    ) # end mainPanel
  ) # end sidebarLayout
  
))
# end shinyUI
  
