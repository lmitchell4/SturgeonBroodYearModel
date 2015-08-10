# These functions were created around Jul/Aug/Sept/Oct 2014 and are used in the
# analyses for mesh selectivity; the mesh being the trammel nets we use for
# sturgeon tagging (6", 7", and 8" sizes); the MeshSelectivity function was
# inspired by the article in MeshSizeSelectivity_Lott1991.pdf

# cleaned up code 10-Jul-2015 & 13-Jul-2015 (J. DuBois) - Note: these functions
# are mostly written to model sturgeon data, attempts to use data other than
# sturgeon with specifically named columns likely will result in errors; in the
# future it would be nice to make these functions more universal to fit other
# data

# Note: copy of this same file (functions_gear_selectivity.R) also can be found
# at C:\Data\..\RSourcedCode (J. DuBois)

DataSubMesh <- function(data, col_date, col_spec,
                        species = c("W", "G"), years) {
  # This function subsets sturgeon data for the purposes of mesh analyses
  
  # Args:
  #   data:     the dataframe containing sturgeon sampling data (or other data)
  #             where sampling gear is involved
  #   col_date: the column which contains the date/time data type
  #   col_spec: the column which contains the species information
  #   species:  the desired species (W = white, G = Green)
  #   years:    the block of years on which to perform subsetting
  
  # Returns:
  #   a dataframe just like data but subsetted accordingly
  
  # Define variables *************************************************
  # get species
  species <- match.arg(arg = toupper(x = species),
                       choices = c("White Sturgeon","Green Sturgeon"),
                       several.ok = TRUE)
  # get date column
  colDate <- eval(expr = substitute(expr = col_date), envir = data,
                  enclos = parent.frame())
  # get species column
  colSpec <- eval(expr = substitute(expr = col_spec), envir = data,
                  enclos = parent.frame())
  
  # locations to remove **********************************************
  # (don't want rescue fish & 1 and 2 are from early years)
  noLoc <- c("Fremont Weir", "Tisdale Bypass", "1", "2")
  
  # mesh size to select on *******************************************
  # hardcoded for now as mesh composition not likely to change
  mesh <- c(6, 7, 8)                                       
  
  # subset data for analysis *****************************************
  new_data <- subset(x = data, 
                     subset = 
                       # gets user-defined years
                       YearFun(x = colDate) %in% years &
                       
                       # removes 2011 rescue data where nets were not used
                       !(Location %in% noLoc) &
                       
                       # gets desired species
                       colSpec %in% species &
                       
                       # selects only mesh size 6,7,8 inch
                       TrammelNetMeshSize %in% mesh       
                     ) #!is.na(LenCat)
  
  # Resetting factors given subsetting *******************************
  rm_species <-
    # gets species to remove
    levels(new_data$Species)[!(levels(new_data$Species) %in% species)]
  
  new_data <- within(data = new_data, expr = {
    Location <- factor(x = Location, exclude = as.factor(noLoc))
    Species <- factor(x = Species, exclude = as.factor(rm_species))
  })
  
  # Function output **************************************************
  #list(Sp = levels(new_data$Species), Loc = levels(new_data$Location))
  new_data
}
# end DataSubMesh

MergeEffortData <- function(data.catch = NULL, data.effort,
                            species = c("W", "G"), years,
                            rm_columns = c("Lat","Lon","Time")) {
  
  # This function merges effort data with catch data and creates one single
  # dataframe
  
  # Args:
  #   data.catch:  NULL (default) will use dfStuAll (all sturgeon data)
  #                otherwise argument should be dataframe that contains
  #                data about sturgeon catch
  #   data.effort: a dataframe that contains data on effort (net sets)
  #                put forth during sturgeon tagging
  #   species:     the desired species (W = white, G = Green)
  #   years:       the block of years on which to perform subsetting
  #   rm_columns:  abbreviation of column names to remove (likely leave default)
  
  # data.catch, species, and years passed to DataSubMesh()
  
  # Returns:
  #   a dataframe containing catch AND effort data
  
  # use dfStuAll (all sturgeon data) if data.catch is null otherwise use
  # argument to data.catch
  data.catch <- if (is.null(x = data.catch)) { dfStuAll } else { data.catch }
  
  # Get catch data (subsetted accordingly)
  data.catch_sub <- DataSubMesh(data = data.catch, col_date = Date,
                                col_spec = Species, species = species,
                                years = years)
  
  # Merge subsetted catch data with effort data and match all with only catch
  # data (the subsetted catch data)
  data.effort <- subset(x = data.effort, subset = YearFun(x = Date) %in% years)
  df_merge <- merge(x = data.effort, y = data.catch_sub, all = TRUE)
  
  # Resetting factors given subsetting (basically resetting location factors)
  uniqLocation <- unique(x = as.character(x = df_merge$Location))
  rm_location <- levels(df_merge$Location)[!(levels(df_merge$Location)
                                             %in% uniqLocation)]
  
  # Reset location factor and change AmtNetOut field to character without /8
  # (without denominator) - done (I think for ease of analyses later on)
  df_merge <- within(data = df_merge, expr = {
    Location <- factor(x = Location, exclude = as.factor(rm_location))
    AmtNetOut <- sub(pattern = "/.", replacement = "",
                     x = as.character(x = AmtNetOut))
  })
  
  # Gives ability to shorten (column-wise) the dataframe (basically allows for
  # removing columns where column name abbreviation is assign to rm_columns
  # parameter)
  rm_columns <- paste(rm_columns, collapse = "|")
  colNames <- grep(pattern = rm_columns, x = colnames(df_merge))
  df_merge <- subset(x = df_merge, select = -colNames)
  
  # Function output
  df_merge
}
# end MergeEffortData

MeshSelectivity <- function(data, len_col, by.break = 5, len_break = NULL,
                            graph = TRUE, ...) {
  
  # This function develops a selectivity curve for sturgeon based on the three 
  # meshes in use since 1990: 6", 7", and 8"; this function will allow for 
  # flexibility in choosing total length versus fork length, all meshes or just 
  # one, location (?), and species; it will also allow flexibility in handling 
  # varying effort between mesh sizes
  
  # Args:
  #   data:      a dataframe of sturgeon data
  #   len_col:   column in dataframe that contains the length data
  #   by.break   value that makes the size of the bins (default 5)
  #   len_break: NULL (default) but can be supplied to change the length bins
  #              for plotting the histogram (density plot)
  #   graph:     TRUE (default) will display selectivity plot
  #   ...:       passed on to col_date, col_species, species, years
  
  # Returns:
  #    a list of selected output - data of proportions and info about
  #    length breaks and bin labels; plus will display moreOutput and
  #    a plot (if desired)
  
  # Subset data
  # ... passed on to col_date, col_species, species, years
  data_sub <- DataSubMesh(data = data, ...)
  
  # Get column of lengths in data_sub
  lenCol <- eval(expr = substitute(expr = len_col), envir = data_sub,
                 enclos = parent.frame())
  
  # Establish length bins and add column to data_sub
  
  # Variables for establishing length bins ***************************
  # ensures smallest length will be included
  startBin <- 10 * floor(x = min(lenCol, na.rm = TRUE) / 10)
  # starts the sequence
  startSeq <- startBin + 25
  # ends the sequence
  endSeq <- 185 #200                                 
  
  # added if/else to allow for different length breaks (e.g., bins for WSTALKEY)
  # - J. DuBois 17-Feb-2014
  if (is.null(len_break)) {
    
    # last bin is endSeq to Inf
    lenBreak <-
      c(startBin, seq(from = startSeq, to = endSeq, by = by.break), Inf)
    # removes the last value (in this case Inf)
    labBreak <- lenBreak[-(length(lenBreak))]
    
  } else {
    
    lenBreak <- len_break
    # removes the last value (in this case Inf)
    labBreak <- lenBreak[-(length(lenBreak))]
    
  }  
  
  # adds a length category field to the dataframe (basically for binning
  # lengths)
  data_sub$LenCat <-
    AddLenCat(col_length = lenCol, breaks = lenBreak, labels = labBreak)
  
  # Adding summaries for data output (column names 'Species', 'Date',
  # 'Location', 'TrammelNetMeshSize' are hard coded, thus this section will
  # throw an error if these column names change within data - something to keep
  # in mind can perhaps deal with this later)
  
  # Possible output to be evaluated later
  moreOutput <-
    substitute(expr =
                 c(
                   tableSpecCount <- with(data = data_sub, expr = {
                     table(YearFun(x = Date), Species, useNA = "ifany")
                   }),
                     
                   tableMeshCatchCount <- with(data = data_sub, expr = {
                     table(YearFun(x = Date), Mesh = TrammelNetMeshSize,
                           useNA = "ifany")
                   }),
                   
                   tableLocationCount <- with(data = data_sub, expr = {
                     table(Location, useNA = "ifany")
                   }),
                   
                   lengthRangeSpecies <- range(lenCol, na.rm = TRUE),
                   
                   dfLenFreqMesh <-
                     ddply(.data = data_sub,
                           .variables = .(Yr = YearFun(x = Date),
                                          Mesh = TrammelNetMeshSize, LenCat),
                           .fun = summarise,
                           Catch = length(x = Species[!is.na(x = Species)]))
                   )
               ) # end substitute
  
  # Remove records where LenCat is na
  data_sub <- subset(x = data_sub, subset = !is.na(LenCat))
  
  # Count all sturgeon LenCat and TrammelNetMeshSize (should colname
  # TrammelNetMeshSize change, I'll need to change below as well)
  dfCount <- ddply(.data = data_sub,
                   .variables = .(LenCat, Mesh = TrammelNetMeshSize),
                   .fun = summarise, .drop = FALSE,
                   nStu = length(Species[!is.na(x = Species)]))
  
  # *** added unbinned 19-Feb-2015 to get a look a the graph **************
  # rather clunky to add it this way, but quick
  dfCountUnBin <- ddply(.data = data_sub,
                        .variables = .(TL_c, Mesh = TrammelNetMeshSize),
                        .fun = summarise, .drop = FALSE,
                   nStu = length(Species[!is.na(x = Species)]))
                   
  dfCountUnBin <- within(data = dfCountUnBin, expr = {
    nStu.new <- ifelse(test = Mesh == "8", yes = nStu / 2, no = nStu)
    Mesh <- factor(x = Mesh) # convert to factor for plotting purposes
  })
  # ***********************************************************************
  
  nMesh <- length(x = unique(x = dfCount$Mesh))
  
  # Divide catch in 8" in half -- see AmtNetOut.R for rationale
  dfCount <- within(data = dfCount, expr = {
    nStu.new <- ifelse(test = Mesh == "8", yes = nStu / 2, no = nStu)
    Mesh <- factor(x = Mesh) # convert to factor for plotting purposes
    # adds lower and upper length range for the purposes of correcting catch
    # data from tagging - may not use below
    #len_low <- as.numeric(x = as.character(x = LenCat))
    #len_up <- c(len_low[-(1:nMesh)] - 1, rep(x = Inf, length.out = nMesh))
  })
  
  # Get proportion of count by LenCat
  dfProp <- ddply(.data = dfCount, .variables = .(LenCat),
                  .fun = mutate, .drop = FALSE,
                  # proportion not corrected for 'effort'
                  pStu = nStu / sum(nStu),
                  # proportion corrected for 'effort'
                  pStu.new = nStu.new / sum(nStu.new)) 
  
  # Plot pStu for each Mesh by LenCat
  labXAxis <- deparse(expr = substitute(expr = len_col))
 
  pltSelectivity <- ggplot(data = dfProp,
                           mapping = aes(x = LenCat, y = pStu.new,
                                         group = Mesh)) + #pStu
    geom_line(mapping = aes(colour = Mesh, linetype = Mesh), size = 1) +
    #geom_point(mapping = aes(fill = Mesh, shape = Mesh),
    #           size = 6, colour = "white") + #, alpha = 1/6
    #scale_fill_manual("",
    #                  values = c("6" = "black", "7" = "red", "8" = "green")) +
    scale_colour_manual("",
                        values = c("6" = "black", "7" = "red", "8" = "green")) +
    scale_linetype_manual("",values = c("6" = "solid", "7" = "twodash",
                                        "8" = "dashed")) +
    #scale_shape_manual("", values = c("6" = 21, "7" = 22, "8" = 24)) +
    #scale_fill_discrete("Mesh") + scale_colour_discrete("Mesh") +  
    labs(
      title = "Sturgeon: Proportion of catch per length category per mesh size",
      y = "Proportion", x = paste("Length category (in ", labXAxis, " cm)")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.50))
          #, size = 10
  
  # forced display of moreOutput
  print(eval(expr = moreOutput))
  
  # Function output
  if (graph) { print(x = pltSelectivity) }
  
  # show some variable in output of list
  list("Output Data" = dfProp,
       "LenRange" = range(data_sub$TL_c, na.rm = T),
       #"Length Range" = lengthRangeSpecies,
       #"Location" = tableLocationCount,
       #"Catch by Mesh" = tableMeshCatchCount,
       #"Catch by Species" = tableSpecCount,
       "Length Break" = lenBreak,
       "Break Labels" = labBreak,
       "UnBinned" = dfCountUnBin # data frame of catch per length not binned
       #levels(cut(x = lenCol, breaks = lenBreak,
       #           right = F, include.lowest = F))
       #ddply(.data = data_sub, .variables = .(Year,Status),
       #      .fun = summarise, Mn = min(TL_c), Mx = max(TL_c))
       )
  #list(lenBreak, labBreak, length(lenBreak)-1) # for testing purposes
  #dfLenFreqMesh
}
# end MeshSelectivity

StuCatchAdjust <- function(data, len_col, years,  by.break,
                           len_break, species = "w") {
  # This function uses MeshSelectivity to adjust catch per mesh from tagging 
  # data; bins are set according to bins in MeshSelectivity, so for finer 
  # adjustment changes in MeshSelectivity need to be made
  
  # Args:
  #   data:      a dataframe of sturgeon catch data with lengths
  #   len_col:   the column in the dataframe that contains the lengths
  #   years:     years on which to subset the dataframe
  #   by.break:  passed to MeshSelectivity()
  #   len_break: passed to MeshSelectivity()
  #   species:   desired species (W (default) = white, G = green)
  
  # Returns:
  #   a summary of catch (as a dataframe) and a force print (display)
  #   of two plots (catch and modeled count vs count)
  
  # Get variables for analyses
  lenCol <- substitute(expr = len_col)
  current_year <- as.numeric(x = format(x = Sys.Date(), format = "%Y"))
  
  # Not elegant, but at least acts as a reminder to change len_col below
  if (!(as.character(x = lenCol) == "TL_c")) {
    stop("len_col argument in 'MeshSelectivity' set at 'TL_c' - 
         improper binning will occur!", call. = FALSE)
  }
  
  # Subset catch data
  dfCatch <- DataSubMesh(data = data, col_date = Date, col_spec = Species,
                         species = species, years = years)
  
  # Call MeshSelectivity (1990 was first year of using three mesh types)
  dfMeshSelect <-
    MeshSelectivity(data = data, len_col = TL_c, by.break = by.break,
                    graph = FALSE, col_date = Date, len_break = len_break,
                    col_spec = Species, species = species, years = 1990:2013)
  #current_year
  
  # Establish variables for setting length category (LenCat) field in dfCatch
  stu_lengths <- eval(expr = lenCol, envir = dfCatch, enclos = parent.frame())
  stu_len_breaks <- dfMeshSelect$`Length Break`
  stu_len_breaks_labs <- dfMeshSelect$`Break Labels`
  
  # Create LenCat field in dfCatch
  dfCatch$LenCat <- AddLenCat(col_length = stu_lengths, breaks = stu_len_breaks,
                              labels = stu_len_breaks_labs)
  
  # Create summary of catch data by year, mesh, and lencat
  # removes where LenCat = NA
  dfCatch <- subset(x = dfCatch, subset = !is.na(x = LenCat))
  dfCatchSummary <- ddply(.data = dfCatch,
                          .variables = .(Yr = YearFun(x = Date),
                                         Mesh = TrammelNetMeshSize, LenCat),
                          .fun = summarise,
                          Count = length(LenCat))
  
  # change Mesh to factor for ease of plotting
  dfCatchSummary$Mesh <- factor(x = dfCatchSummary$Mesh)
  
  # With the two dataframes below where levels(Mesh) and levels(LenCat) match, 
  # get pStu.new from dfMeshSelect$`Output Data`
  
  #dfMeshSelect$`Output Data`
  #dfCatchSummary
  
  # Get row index from dfMeshSelect$`Output Data` where mesh and lencat match
  # with dfCatchSummary
  indexMatch <-
    match(x = paste(dfCatchSummary$Mesh, dfCatchSummary$LenCat, sep = "_"),
          table = paste(dfMeshSelect$`Output Data`$Mesh,
                        dfMeshSelect$`Output Data`$LenCat, sep = "_")
          )
  
  # Add proportion to dfCatchSummary and adjust count (catch)
  dfCatchSummary <- within(data = dfCatchSummary, expr = {
    PropAdj <- dfMeshSelect$`Output Data`$pStu.new[indexMatch]
    CountAdj <- Count / PropAdj
  })
  
  # Plot pStu for each Mesh by LenCat
  labXAxis <- deparse(expr = substitute(expr = len_col))
  
  # sets up y-position for plotting text on graph
  yPos <- with(data = dfCatchSummary, expr = max(CountAdj, na.rm = TRUE)) 
  
  # create dataframe of stats to be displayed for information
  dfStats <- 
    ddply(.data = dfCatchSummary, .variables = .(Yr, Mesh),
          .fun = here(f = summarise),
          Rsq = round(x = summary(lm(formula = CountAdj ~ Count))$r.squared,
                      digits = 2),
          Slp = round(x = coef(lm(formula = CountAdj ~ Count))[2], digits = 2), 
          nAdj = round(x = sum(CountAdj), digits = 0),
          nNot = sum(Count),
          # x position of text on plot levels(LenCat)[max(LenCat)]
          # 26 = max level
          xVal = 26,
          # y position of text on plot
          yVal = floor(x = yPos / 10) * 10
          )
  
  # plotting
  set_alpha <- 5/5 # alpha variable for overplotting (if desired)
  pltCatch <- ggplot(data = dfCatchSummary, mapping = aes(x = LenCat)) +
    geom_bar(mapping = aes(y = CountAdj, fill = "black"), stat = "identity") +
    geom_bar(mapping = aes(y = Count, fill = "red"), stat = "identity",
             alpha = set_alpha, width = 0.65) +
    # can replace Mesh with . for annual LF distributions
    facet_grid(facets = Yr ~ Mesh) + 
    scale_fill_manual(name = "Catch", guide = "legend", labels = c("Adj","Not"),
                      values = c("black"="black","red"="red")) +
    geom_text(mapping = aes(x = xVal, y = yVal,
                            label = paste0("Not=",nNot,"\n","Adj=",nAdj)),
              data = dfStats, size = 3, hjust = 1, vjust = 1) +
    labs(title = "Sturgeon: Catch and adjusted catch per mesh size",
         y = "Catch", x = paste("Length category (in ", labXAxis, " cm)")) +
    theme(axis.text.x = element_text(angle = 90, size = 12,
                                     hjust = 1, vjust = 0.35),
          legend.position = "bottom")
    # hjust = 1; vjust = 0.15 (currently)
  
  pltScatter <- ggplot(data = dfCatchSummary,
                       mapping = aes(x = Count, y = CountAdj)) +
    stat_smooth(method=lm, se = TRUE) +
    geom_point() +
    facet_grid(facets = Yr ~ Mesh) +
    geom_text(mapping = aes(x = 2, y = 145, label = paste0("R^2: ", Rsq)),
              data = dfStats, parse = TRUE, size = 3) +
    geom_text(mapping = aes(x = 10, y = 145, label = paste0("Slope: ", Slp)),
              data = dfStats, parse = FALSE, size = 3)
  
  # Summary of catch just by year
  dfStatsTotal <- 
    ddply(.data = dfCatchSummary, .variables = .(Yr), .fun = summarise,
          nStu = sum(Count),
          nAdj = round(x = sum(CountAdj), digits = 0)
          )
  
  # Summary of catch & adjusted catch by LenCat for each year
  dfCatchByLenCatSum <- 
    ddply(.data = dfCatchSummary, .variables = .(Yr, LenCat), .fun = summarise,
          # straight catch
          n = sum(Count),
          # catch adjusted by proportions from mesh selectivity
          n.adj = sum(CountAdj)
          )
  
  # Function output
  print(x = pltCatch)
  print(x = pltScatter)
  dfStatsTotal
  dfCatchSummary
  #head(dfCatchSummary[,c(1,2,3,4,6,5)])
  #levels(dfCatchSummary$LenCat)
  #dfCatchByLenCatSum
  #head(dfCatchSummary)
  #dfCatchSummary
}
# end StuCatchAdjust

MeshCatchByYear <- function(data = dfStuGearData_sub, year) {
  # This function calculates catch per mesh size for the appropriate year; data
  # is subsetted already to give only years 1990:2013 and to remove any records
  # where length is NA
  # ***************************************************************************
  
  # Args:
  #   data: previously subsetted dataframe where any records for length = NA
  #         have been removed
  #   year: year on which to subset data (vector with length = 1)
  
  # Returns:
  #   a dataframe of count of each length for each mesh size and where mesh
  #   size was not given (0 or NA)
  
  # Subset data
  dfYear <- subset(x = data, subset = Year %in% year)
  
  # Calculate catch per mesh size (6", 7", and 8" panels)
  dfCountMesh <- ddply(.data = dfYear, .variables = .(TL_c), .fun = summarise,
                       # added MeshNA 26-Nov-2014 to get an idea of data not
                       # being used in analyses (note some mesh data entry as 0
                       # synonymous with NA) - J. DuBois
                       MeshNA = length(x = Species[TrammelNetMeshSize %in%
                                                     c(0, NA)]),
                       Mesh6 = length(x = Species[TrammelNetMeshSize %in% 6]),
                       Mesh7 = length(x = Species[TrammelNetMeshSize %in% 7]),
                       Mesh8 = length(x = Species[TrammelNetMeshSize %in% 8])
                       )
  
  # Function output
  dfCountMesh
}
# end MeshCatchByYear

ExpandDataAppFun <- function(dtfrm, col_data, col_reps, FUN = mean,
                             transform = FALSE) {
  # This function applies a function (e.g., mean) to a dataframe column 
  # (col_data) that has been expanded according to data (e.g., catch) in another
  # column or columns; the data on which the function is being employed can be 
  # log (natural) transformed if desired; function output is the applied 
  # functions output per each of the repeated columns (cols_rep); best currently
  # applied to a catch table by length where mean length is desired on expanded 
  # catch; designed sort of like what is described on p. 189 of The R Book 
  # ************************************************
  
  # Args:
  #   dtfrm:     a dataframe on which to apply the function (FUN)
  #   col_data:  column on which to apply fun and transform if desired
  #   col_reps:  column that contains the (count) data used to expand
  #              each length value
  #   FUN:       the function to be applied to the data (default: mean)
  #   transform: natural log transform the data if desired (default: FALSE)
  
  # Return:
  #   results of the function applied to the selected data (column)
  
  # get function according to input
  fun <- match.fun(FUN = FUN)
  
  # get data column on which to apply fun and transform if desired
  dat <- if (transform) {
    log(x = dtfrm[, col_data]) # natural log (base = exp(1))
    } else {
      dtfrm[, col_data]
    }

  # created if/else to properly apply function if only col_reps was only 1 
  # column
  if (length(col_reps) > 1) {
    
    # set dataframe columns as a list
    cols_to_rep <- as.list(dtfrm[, col_reps])
    
    # apply selected function (fun) over list
    res <- sapply(X = cols_to_rep, FUN = function(cols) {
      fun(rep(x = dat, times = cols))
    }, simplify = TRUE, USE.NAMES = TRUE)
    
  } else {
    
    # single column selected
    cols_to_rep <- dtfrm[, col_reps]
    
    # apply selected function (fun) to data
    res <- fun(rep(x = dat, times = cols_to_rep))
    
    # set name for output as single column name
    names(res) <- colnames(x = dtfrm[col_reps])
    
  }
  
  # function output
  res
}
# end ExpandDataAppFun

MeshCatch <- function(stu_data = dfStuAll, year_end, min_len = 0,
                      even_mesh_only = FALSE, plot_den = FALSE) {
  # This function calculates the catch per mesh for subsetted sturgeon data; 
  # data subsetted for species (white), years (1990 to year_end), and minimum
  # length (total length)
  # ************************************************
  
  # Args:
  #   stu_data:       a dataframe of sturgeon data with length column   
  #   year_end:       end year on which to subset data (will begin and 1990)
  #   min_len:        minimum length on which to subset and then (eventually)
  #                   model (defaults to include all; set to 0)
  #   even_mesh_only: option to only use amount of net out of 4/8 or 8/8
  #                   (default: FALSE)
  #   plot_den:       option to display a density plot for lengths
  #                   (default: FALSE)
  
  # Returns:
  #   a list where each element is a dataframe (of counts and of the subsetted
  #   data)
  
  # stu_data subset
  dfStu <- DataSubMesh(data = stu_data, col_date = Date, col_spec = Species,
                       species = "w", years = 1990:year_end)
  
  # further subsetting: removes records where TL_c is NA and inlcudes only those
  # meeting minimum length (min_len) requirement
  dfStu_sub <- subset(x = dfStu, subset = !is.na(x = TL_c) & TL_c >= min_len)
  
  # further subsetting: option to only use amount of net out of 4/8 or 8/8 (thus
  # using only data where ratio of 6:7:8 mesh was 1:1:2)
  dfStu_sub <- if (even_mesh_only) {
    
    # get dates from dfStuEffort where amount net out = 4/8 or 8/8
    dates <- with(data = dfStuEffort, expr = {
      unique(x = Date[as.character(x = dfStuEffort$AmtNetOut)
                      %in% c("4/8", "8/8")])
    })
    
    # subset on dates
    dfStu_sub[dfStu_sub$Date %in% dates,]
    
    } else {
      dfStu_sub
    }
  
  # add year field (as numeric) to dataframe (dfStu_sub)
  dfStu_sub$Year <- YearFun(x = dfStu_sub$Date, output = "asNum")
  
  # create dataframe of catch by length for each mesh size
  dfCountMesh <-
    ddply(.data = dfStu_sub, .variables = .(TL_c), .fun = summarise,
          Mesh6 = length(x = Species[TrammelNetMeshSize %in% 6]),
          Mesh7 = length(x = Species[TrammelNetMeshSize %in% 7]),
          Mesh8 = length(x = Species[TrammelNetMeshSize %in% 8])
          )
  
  # option for density plot on raw data
  if (plot_den) {
    plt <- ggplot(data = dfStu_sub,
                  mapping = aes(x = TL_c, fill = factor(TrammelNetMeshSize))) +
      geom_density(alpha = 2/5) +
      labs(fill = "Mesh", x = "Length (cm TL)")
    # forced display of plot
    print(plt)
  }
  
  # function output
  list(CatchData = dfCountMesh, Data = dfStu_sub)
}
# end MeshCatch

RelRetStartVals <- function(catch_data, mesh_size) {
  # This function gets the potential starting values for parameters used in 
  # Millar code (selectivity curves algorithms)
  
  # Args:
  #   catch_data: dataframe containing column of sturgeon catch
  #   mesh_size: the mesh size (in cm)  for each net panel
  
  # Returns:
  #    a dataframe of summary data and a list containing the range for each
  #    column in the summary data
  
  # NOTE: may want to consider adding some kind of validation for centimeters
  # versus inches for the argument mesh_size
  
  # ************************************************
  
  # get data (not sure why I did this, but will keep it for now - J. DuBois 
  # 13-Jul-2015)
  mesh_catch_data <- catch_data
  
  # variables for this function
  ms <- mesh_size # converted to cm; use same units for fish length and mesh
  relative_mesh_size <- ms / min(ms)
  
  # variables for ExpandDataAppFun function
  col_data = "TL_c"
  col_reps = 2L:4L
  
  # list used for ExpandDataAppFun funtion (transform = TRUE means data will be
  # natural log transformed prior to analyses)
  lstData <- list(mean_length    = c(FUN = "mean", transform = FALSE),
                  sd_length      = c(FUN = "sd", transform = FALSE),
                  mean_length_lg = c(FUN = "mean", transform = TRUE),
                  sd_length_lg   = c(FUN = "sd", transform = TRUE)
                  )
  
  # apply ExpandDataAppFun function to lstData
  lstDataSummary <- sapply(X = names(lstData), FUN = function(x) {
    ExpandDataAppFun(dtfrm = mesh_catch_data, col_data = col_data,
                     col_reps = col_reps,
                     FUN = lstData[[x]][1], transform = lstData[[x]][2])
  }, simplify = TRUE, USE.NAMES = TRUE)
  
  # convert list to dataframe
  dfDataSummary <- as.data.frame(x = lstDataSummary)
  
  # add columns to dataframe
  dfDataSummary$k_guess <- dfDataSummary$mean_length / relative_mesh_size
  dfDataSummary$k_guess_lg <- dfDataSummary$mean_length_lg / relative_mesh_size
  
  # get ranges of values for each column in dfDataSummary
  lstDataRanges <- sapply(X = dfDataSummary, FUN = range,
                          simplify = TRUE, USE.NAMES = TRUE)
  
  # function output
  list(Data = dfDataSummary, Ranges = lstDataRanges)
}
# end RelRetStartVals

ToList <- function(dtfrm = mtrx_seq, fifth_p = fp_start_val, rn) {
  # This function takes the dtfrm input and creates a list on each row of the 
  # dtfrm; row selected per argument 'rn; creating a list aids in running the
  # Millar model fit functions using multiple parameters
  
  # Args:
  #   dtfrm:   the dataframe of parameters for the Millar model algorithms
  #   fifth_p: the fifth parameter for binorm and bilognorm curve fits 
  #   rn:      the rownumer in dtfrm of the desired set of parameters
  
  # Returns:
  #   a list of parameters for the Millar model curve fits (algorithms)
  
  # get number of rows in dtfrm
  num_rows <- nrow(x = dtfrm)
  
  # stop function if 'rn' value is not a row number of the dtfrm
  if (rn < 1 || rn > num_rows) {
    stop(paste0("Row ", rn, " does not exist in dataframe: ",
                substitute(dtfrm)), call. = FALSE)
  }
  
  # combine the selected row (by row number) into a list
  res <-  list(
    norm.loc = c(dtfrm[rn, "avg1"], dtfrm[rn, "sd1"]),
    norm.sca = c(dtfrm[rn, "avg1"], dtfrm[rn, "sd1"]),
    lognorm = c(dtfrm[rn, "avglog1"], dtfrm[rn, "sdlog1"]),
    binorm.sca = c(dtfrm[rn, "avg1"], dtfrm[rn, "sd1"],
                   dtfrm[rn, "avg2"], dtfrm[rn, "sd2"], fifth_p),
    bilognorm = c(dtfrm[rn, "avglog1"], dtfrm[rn, "sdlog1"],
                  dtfrm[rn, "avglog2"], dtfrm[rn, "sdlog2"], fifth_p)
    )
  
  # function output
  res # res is list
}
# end ToList

RelRetModelFit <- function(catch_data, val_list, plotlens = NULL,
                           mesh_size = mesh_size_cm, win = FALSE,
                           plot_n_data = FALSE) {
  # This function applies Millar code (selectivity curve models) to catch per
  # mesh data
  
  # Args:
  #   catch_data:  a dataframe of sturgeon data with column of catch
  #   val_list:    value list of which to apply Millar net fit models
  #   plotlens:    desired length range on which to plot Millar models using
  #                PlotCurves function
  #   mesh_size:   the mesh size (in cm) of each net panel
  #   win:         gives option for displaying plots in separate plot window;
  #                these plots display the output of fitting the Millar
  #                model types
  #   plot_n_data: (added 28-Jul-2015) adds option to show deviance and gear
  #                selectivity plots and data; wanted this option primarily
  #                for use with brood model and not wanting to see Millar
  #                ouput with brood model output
  
  # Returns:
  #    (invisible) list of curve data, value lists (val_list), and summary data
  #    plus a plot of the model deviance and model curves
  # ************************************************
  
  # establish fishing effort (hard coded for now, but could as an argument 
  # later)
  mesh_rel_power <- c(0.25, 0.25, 0.50) # (theoretically) should sum to 1
  
  # establish fifth parameter for binorm.sca and bilognorm curve fits
  #fp_start_val <- fifth_para_start_val # default of 0.65 chosen arbitrarily
  #dblFifthParam <- log(fp_start_val / (1 - fp_start_val))
  
  # curve list on which to apply Millar net fit models
  curve_list <- list(norm.loc = "norm.loc", norm.sca = "norm.sca",
                     lognorm = "lognorm", binorm.sca = "binorm.sca",
                     bilognorm = "bilognorm"
                     )
  
  # create value list of which to apply Millar net fit models
  value_list <- val_list
  #list(
  #norm.loc = avg_sd_len,
  #norm.sca = avg_sd_len,
  #lognorm = avg_sd_log_len,
  #binorm.sca = c(avg_sd_len, avg_sd_len2, dblFifthParam),
  #bilognorm = c(avg_sd_log_len, avg_sd_log_len2, dblFifthParam)
  #) 
  
  # for output format (removing for now, J. DuBois 15-Jul-2015)
  #cat("\n********* output separator: below run with different starting values ",
      #deparse(substitute(expr = val_list)),
  #    "******************\n", sep = "") 
  
  # for output format
  #cat("\n****************** Model fit parameters & deviance ****************\n")
  
  # apply Millar 'NetFit' funtion to data for each model type (in curve_list)
  fit_list <- mapply(FUN = NetFit, rtype = curve_list, x0 = value_list,
                     MoreArgs = list(Data = catch_data, Meshsize = mesh_size,
                                     rel.power = mesh_rel_power,
                                     display_param = plot_n_data),
                     SIMPLIFY = FALSE
                     )
  
  #cat("*********************************************************************\n")
  
  # establish estimates variable with Millar 'Estimates' function
  estimates <- sapply(X = fit_list, FUN = Estimates,
                      simplify = TRUE, USE.NAMES = TRUE)
  # Shows table with Deviance and degrees of freedom (d.o.f); ideally want low 
  # deviance and deviance to be around the d.o.f value (see bottom of page 102 
  # Millar); also displays deviance plots
  
  # set up plot grid for deviance plots and relative retention plots
  if (win) {
    # option to turn this feature off
    windows(width = 8, height = 12) # plot window for monitor view
  }
  
  par(mfcol = c(5, 2), mar = c(4.1, 4.1, 1, 2))
  
  # apply Millar 'Summary' function to model fit data
  summary_data <- data.frame(variable = c("null.l", "model.l", "full.l",
                                          "Deviance", "Pearson.chisq", "d.o.f"),
                             sapply(X = fit_list, FUN = Summary,
                                    show_plots = plot_n_data,
                                    simplify = TRUE, USE.NAMES = TRUE)
                             )
  
  # apply Millar 'PlotCurves' function to data; displays relative retention 
  # plots and data behind plots
  curve_data <- sapply(X = fit_list, FUN = PlotCurves, plotlens = plotlens,
                       show_plots = plot_n_data,
                       simplify = FALSE, USE.NAMES = TRUE)
  
  # for output format
  #cat("\n*************** Model fit summary data ***************\n\n")
  #print(curve_data)
  
  # function output
  # display estimates and summary data
  #print(list(Estimates = estimates, SummaryTable = summary_data))
  # data used for further testing when needed (thus invisible())
  #invisible(x = c(curve_data, value_list))
  invisible(list(CurveData = curve_data, ValueList = value_list,
                 Estimates = estimates, SummaryTable = summary_data)
            ) # changed output 09-Oct-2014 @ 1535
}
# end RelRetModelFit

RelRetAdjCatch <- function(dtfrm, model, env = model_fit) {
  # This function adjusts catch from sturgeon tagging using relative retention 
  # values from mesh selectivity analyses (using Millar code); idea is to plot 
  # (boxplot) all length values (from adjusted catch) as function of Year, plot 
  # should show areas of good recruitment (as indicated by a lower median or 
  # 25th percentile)
  
  # Args:
  #   dtfrm: a dataframe of sturgeon catch data
  #   model: a chosen (Millar) model on which to model (adjust) catch
  #   env:   the environment in which to evaluate the model
  
  # Returns:
  #   a list with columns of catch and adjusted catch totals
  # ************************************************
  
  # get unique years from data
  # Year should be a field in dtfrm, if not error
  Years <- as.list(x = unique(x = dtfrm$Year))
  
  #unlist(x = Years, use.names = FALSE)
  names(Years) <- unique(x = dtfrm$Year) # list named with years
  
  # Get data from each year in Years
  lst_catch <- sapply(X = Years, FUN = MeshCatchByYear, data = dtfrm,
                      simplify = FALSE, USE.NAMES = TRUE)
  
  # Data from model fit; envir = model_fit is global variable, if this changes 
  # must change here too (Note: removed eval() statement in favor of 
  # env[["CurveData"]][model] to allow for easier changing of model type say 
  # from "bilornorm" to "binorm". "CurveData" hardcoded as it comes from
  # invisible output of RelRetModelFit - J. DuBois 13-Jul-2015)
  model_fit_data <- env[["CurveData"]][model]
  #eval(expr = substitute(expr = model), envir = env)
  
  # convert to data frame for analyses
  curve_data <- as.data.frame(x = model_fit_data) 
  
  # Gets row number from dfCurveData where length (TL_c) in lstCatchData matches
  # length (TL_c) in dfCurveData
  length_match <- sapply(X = Years, FUN = function(y) {
    y <- as.character(y)
    match(x = lst_catch[[y]][["TL_c"]], table = curve_data[,1])
  }, simplify = TRUE, USE.NAMES = TRUE)
  
  # Creates list (by year) with length, length category, catch, and adjusted
  # catch data
  lst_adj_catch <- mapply(FUN = function(x, y) {
    y <-  as.character(y)
    cbind(Len = lst_catch[[y]][, "TL_c"],
          # gets length data for each year from lst_catch adds length category
          # (as ordered factor) for next summary steps; NOTE LenCat add category
          # sub, leg, ovr based on <= 116 sub, 117-168 leg, >= 169 ovr Should
          # category change to based on status at tagging?? (01-Dec-2014)
          
          LenCat = ordered(x = LenCat(x = lst_catch[[y]][, "TL_c"]),
                           list("sub" = "sub", "leg" = "leg", "ovr" = "ovr")),
          # removed (for now) adding original catch data +1 added for column 
          # adjustment (J. DuBois, 26-Nov-2014); NOTE not outputting here the 
          # column that holds MeshNA (in this context MeshNA will all be 0
          # because of the way data is subsetting in DataSubMesh function)
          lst_catch[[y]][c(2, 3, 4) + 1],
          
          # gets catch data for each year from lstCatchData; divides catch for
          # each year by matched relative retention value (adjusted catch)
          adj = lst_catch[[y]][c(2, 3, 4) + 1] / curve_data[x, -1])
  }, x = length_match, y = Years, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  # Creates total columns for catch and adjusted catch
  lst_totals <- sapply(X = Years, FUN = function(y) {
    # Create function to add total (for each year and by length) catch and
    # adjusted catch
    
    y <- as.character(y) # set year to character for subsetting below
    
    # assign data to variable
    catch_data <- lst_adj_catch[[y]]
    
    # combine data with catch and adj catch total columns
    cbind(
      # data
      catch_data,
      # catch total
      Catch = catch_data[, "Mesh6"] + catch_data[, "Mesh7"] +
        catch_data[, "Mesh8"],
      # adjusted catch total
      AdjCatch = catch_data[, "adj.Mesh6"] + catch_data[, "adj.Mesh7"] +
        catch_data[, "adj.Mesh8"]
      )
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # Function output
  lst_totals
}
# end RelRetAdjCatch

PlotAdjCatch <- function(lst = lstAdjCatch, catch_type = "AdjCatch",
                         sum_data = FALSE) {
  # This function takes the data from lstAdjCatch and expands it (using 'rep' 
  # function) to in turn plot a boxplot of length as a function of year (tagging
  # year); thinking (according to K. Newman [USFWS]) is that year with lower 
  # median or 25th percentile would indicate a year of good recruitment
  
  # Args:
  #   lst:        a list of data with modeled catch and raw catch for each
  #               tagging year
  #   catch_type: which catch to plot - modeled (AdjCatch) or raw (Catch)
  #   sum_data:   option to include data summary (by year) with plot
  
  # Returns:
  #   a box plot where x-axis is year (tagging year) and y-axis is length
  #   (either fork or total) to show variation in median lenght and range
  #   of length over time series
  
  # see Years as list for further testing below; names for the incoming list (in
  # 'lst') should be the representative tagging year
  Years <- names(lst)
  
  # using the 'rep' function, expand the catch number for each length; this is 
  # to develop a dataframe of Year and length where every length is accounted 
  # for (e.g., 2012 had 8 fish at 85 cm TL, should produce 8 rows of 2012 [Year 
  # column] and 85 [length column])
  lst_expanded_length <- sapply(X = Years, FUN = function(y){
    # AdjCatch or Catch
    rep(x = lst[[y]][, "Len"], times = lst[[y]][, catch_type])
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # convert list above to dataframe
  dfExpandedLength <- ldply(.data = lst_expanded_length, .fun = cbind,
                            .id = "Year")
  
  colnames(dfExpandedLength)[2] <- c("Len")
  
  # convert Year field in dataframe from factor to numeric
  dfExpandedLength$Year <-
    as.numeric(x = as.character(x = dfExpandedLength$Year))
  
  # calculated for each year the catch
  data_sum <- ddply(.data = dfExpandedLength, .variables = .(Year),
                    .fun = summarise,
                    Catch = length(x = Len[!is.na(x = Len)])
                    )
  
  # set up minor breaks for minor break lines in boxplot below
  minor_brks <- seq(from = as.numeric(Years[1]),
                    to = as.numeric(Years[length(Years)]), by = 1)
  
  # create ggplot2 boxplot (_hline values denote current legal-sized slot limit 
  # of 117-168 cm TL; actually slot is now in FL (102-152), but 117-168 is the 
  # total length equivalent)
  bx_plot <- ggplot(data = dfExpandedLength,
                    mapping = aes(x = Year, y = Len, group = Year)) +
    geom_hline(yintercept = c(117, 168), colour = "red", size = 0.5,
               linetype = 2) +
    geom_boxplot(outlier.size = 1) +
    stat_summary(fun.y = mean, geom = "point", colour = "blue", size = 2) +
    scale_x_continuous(minor_breaks = minor_brks, expand = c(0, 0.3)) +
    annotate(geom = "text", x = 1990, y = 275, label = "blue dot = mean",
             size = 3, hjust = 0) +
    annotate(geom = "text", x = 1990, y = 265,
             label = "red line = lower/upper slot", size = 3, hjust = 0) +
    labs(title = catch_type, y = "Length (cm TL)") +
    theme(plot.title = element_text(size = 10))
  
  # function output with option to print plot and data or just plot
  if (sum_data) {
    print(data_sum)
    print(bx_plot)
  } else {
    print(bx_plot)
  }
}
# end PlotAdjCatch

GetGearSelectModeledCatch <- function(min_len = 85, plot_den = FALSE,
                                      millar_model = "bilognorm",
                                      param_row_select = 3, win = FALSE,
                                      max_year = NULL, plot_n_data = FALSE) {
  # function used to create list of catch and adjust catch by year; defaults set
  # for using 70 cm TL as minimum cutoff length; J. DuBois (15-Dec-2014)
  
  # Args:
  #   min_len:     minimum* total length on which apply the Millar model fits
  #                default 85 (cm TL)
  #   plot_den:    option for which to display density plot of mesh catch
  #                default = FALSE
  #   plot_n_data: option to show (or not) Millar plots and data output
  #                added 28-Jul-2015 (J. DuBois)
  
  # *disclaimer: The default is set so as not to include fish where catch at
  # length was exceptionally low (say < 5); including small fish where catch is
  # low tends to cause gear selectivity to be low (for that length) which in
  # turn sends through the roof the modeled catch for that length
  
  # Returns:
  #   nothing; but assigns 'lstAdjCatch' to the global environment;
  #   'lstAdjCatch' contains raw catch and modeled catch by length for each
  #   tagging year
  
  # if list lstAdjCatch exists stop and show warning
  if (exists("lstAdjCatch")) {
    stop("lstAdjCatch already exists. Please check output.", call. = FALSE)
  }
  
  # use current year if max_year is NULL (default)
  if (is.null(max_year)) {
    max_year <- YearFun(Sys.Date(), output = "asNum")
  }
  
  # Step 1 - Get catch data and starting values -----------------------------
  
  # Get data (catch per mesh); data for white sturgeon from 1990 to 2013
  
  # MeshCatch function output also includes subsetted raw data (output is a
  # list: "CatchData", "Data") argument defaults: stu_data = dfStuAll, min_len =
  # 0, even_mesh_only = FALSE, plot_den = FALSE
  dfMeshCatch <- MeshCatch(year_end = max_year, min_len = min_len,
                           plot_den = plot_den)
  
  # Get mesh size data (convert to cm; use same units for length and mesh)
  mesh_size_cm <- c(6, 7, 8) * 2.54
  
  # Get potential starting values for Millar selectivity code (model fits)
  start_values <- RelRetStartVals(catch_data = dfMeshCatch$CatchData,
                                  mesh_size = mesh_size_cm)
  
  # Step 2 - create  matrix of starting values for easier testing -----------
  
  # set starting values for mtrx_seq
  start_range <- start_values$Ranges[, c(5, 2, 6, 4)]
  
  # apply function to selected columns of mean and sd ranges
  mtrx_seq <- apply(X = start_range, MARGIN = 2, FUN = function(x) {
    # function creates a lookup table for variable used in floor/ceiling 
    # functions; different divisors/multipliers needed to get correct range when
    # using floor/ceiling functions
    
    # create look-up data frame
    nchar_lkp <- data.frame(n_char = c(4, 3, 2, 1), val = c(10, 10, 1, 0.1))
    
    # create divisor/multiplier variable for low and high range values; probably
    # could use the same variable for both low and high but made a separate one
    # for each just to be safe and clear
    div_low <- nchar_lkp$val[match(x = nchar((x[1] * 10) %/% 1),
                                   table = nchar_lkp$n_char)]
    div_high <- nchar_lkp$val[match(x = nchar((x[2] * 10) %/% 1),
                                    table = nchar_lkp$n_char)]
    
    # create low and high values on which to sequence
    low <- floor(x = x[1] / div_low) * div_low
    high <- ceiling(x = x[2] / div_high) * div_high
    
    # sequence from low to high for each column; length.out of 5 arbitrary but
    # felt 5 was enough values on which to test
    res <- seq(from = low, to = high, length.out = 5)
    
    # function output
    res # res is a matrix
  })
  
  # rename maxtrix columns
  colnames(mtrx_seq) <- c("avg1", "sd1", "avglog1", "sdlog1")
  
  # add additional columns to dataframe (these columns contain the 3rd and 4th 
  # paramters for the models that require them); column multipliers selected 
  # arbitrarily and can be modified accordingly; mtrx_seq converts to dataframe
  # at this step
  mtrx_seq <- within(data = data.frame(mtrx_seq), expr = {
    sdlog2 <- sdlog1
    avglog2 <- avglog1 * 1.2
    sd2 <- sd1 * 1.1
    avg2 <- avg1 * 1.75
  })
  
  # get fifth parameter for binorm and bilognorm curve fits
  seed_val <- 0.65 # default of 0.65 chosen arbitrarily
  fp_start_val <- log(seed_val / (1 - seed_val))
  
  # Step 3 - Get model fit (Millar code application) ------------------------
  
  # lengths on which to produce relative retention values
  plot_lens <- seq(from = 40, to = 300, by = 1)
  
  model_fit <- RelRetModelFit(catch_data = dfMeshCatch$CatchData,
                              val_list = ToList(dtfrm = mtrx_seq,
                                                fifth_p = fp_start_val,
                                                rn = param_row_select),
                              mesh_size = mesh_size_cm,
                              plotlens = plot_lens, win = win,
                              plot_n_data = plot_n_data)
  
  # Step 4 - Use chosen model fit to adjust catch ---------------------------
  
  # assign a variable to 'RelRetAdjCatch'; output is a list and list will be
  # used to plot time series of boxplot (total length ~ year); change model
  # argument accordingly
  lstAdjCatch <- RelRetAdjCatch(dtfrm = dfMeshCatch$Data,
                                model = millar_model,
                                env = model_fit)
  
  # function output
  #assign(x = "lstAdjCatch", value = lstAdjCatch, envir = .GlobalEnv)
  
  # Note: for now (14-Jul-2015) decided to comment out the assign (lstAdjCatch) 
  # command; rather I'll add to list below the data as (converted to) dataframe 
  # this will make for easier changing of data say when changing some of this
  # function's arguments (J. DuBois)
  
  # convert list to dataframe and then convert Year from factor to character
  df_data <- plyr::ldply(.data = lstAdjCatch, .id = "Year")
  df_data$Year <- as.character(df_data$Year)
  
  # still toying with the output (13-Jul-2015, J. DuBois)
  list(
    # convert lstAdjCatch to dataframe for output
    Data = df_data,
    Model = millar_model,
    Summary= model_fit$SummaryTable,
    #Estimates = model_fit$Estimates,
    MatrixData = mtrx_seq,
    Years = names(lstAdjCatch)
    # an option for Years: noquote(paste0(names(lstAdjCatch), collapse = ", "))
  )
}
# end GetGearSelectModeledCatch

GetStuTagGearSelData <- function(catch_col = c("Catch", "AdjCatch"),
                                 min_len = 85, millar_model = "bilognorm",
                                 param_row_select = 3, l_mit_rows = FALSE,
                                 plot_n_data = FALSE) {
  # This function summarizes count (catch) by year and fork length (FL)
  
  # Args:
  #   catch_col:        the column in the dataframe that contains the catch
  #                     (count) info; can be either raw (Catch) or modeled
  #                     (AdjCatch)
  #   min_len:          the minimum length on which to perform the gear
  #                     selectivity models
  #   millar_model:     the desired Millar model (list here:)
  #   param_row_select: the row in the list of parameters for the Millar model;
  #                     default = 3 (for more info ??????????)
  #   l_mit_rows:       return dataframe with rownames like L. Mitchell who
  #                     uses YYYY.FL notation (default: FALSE, use 1:end)
  
  # Returns:
  #   a dataframe of sum of catch (count) by year by fork length (FL or FL_cm)
  
  # parameters with defaults for GetGearSelectModeledCatch() function
  # min_len = 85                    will keep in automation
  # plot_den = FALSE                will use default
  # millar_model = "bilognorm"      will keep in automation
  # param_row_select = 3            will keep in automation
  # win = FALSE                     for now (14-Jul-2015) will set to TRUE
  # max_year = NULL                 will use default (current year)
  
  # match argument with a selections
  catch_col_select <- match.arg(arg = catch_col,
                                choices = c("Catch", "AdjCatch"),
                                several.ok = FALSE)
  
  # get data from gear selectivity function (uses Millar functions)
  df_catch_data <-
    GetGearSelectModeledCatch(min_len = min_len, millar_model = millar_model,
                              param_row_select = param_row_select,
                              win = FALSE, plot_n_data = plot_n_data)$Data
  
  # add column of fork length (since Len column is total length)
  # algorithm for conversion from biological opinion
  df_catch_data$FL_cm <- round((0.9036 * df_catch_data$Len) - 1.2162,
                               digits = 0)
  
  # select only Year, fork length, and catch or modeled (adj) catch columns
  df_catch_data_sub <- df_catch_data[, c("Year", "FL_cm", catch_col_select)]
  
  # change column names for consistency with Mitchell Model and for aggregate
  # function below
  colnames(df_catch_data_sub) <- c("Year", "ForkLength", "Count")
  
  # sum catch by year - slightly different than L. Mitchell: I sum first then 
  # round, she rounded first then summed; convert Count to integer for
  # consistency with Card data
  df_catch_sum <- aggregate(formula = Count ~ Year + ForkLength,
                            data = df_catch_data_sub,
                            FUN = function(x) {
                              #as.integer(round(sum(x, na.rm = TRUE),
                              #                 digits = 0))
                              # round then sum (as L. Mitchell)
                              as.integer(
                                sum(round(x, digits = 0), na.rm = TRUE)
                              )
                            })
  
  # convert year to numeric for use in L.Mitchell Model
  df_catch_sum$Year <- as.numeric(df_catch_sum$Year)
  
  # re-order (Year, FL) for convenience
  df_catch_sum <- df_catch_sum[order(df_catch_sum$Year,
                                     df_catch_sum$ForkLength), ]
  
  # renumber for convenience (1 thru end) or like L. Mitchell
  rownames(df_catch_sum) <- if (l_mit_rows) {
    paste(df_catch_sum$Year, df_catch_sum$ForkLength, sep = ".")
  } else {
    1:nrow(df_catch_sum)
  }
  
  # function output (summary dataframe)
  df_catch_sum
}
# end GetStuTagGearSelData

