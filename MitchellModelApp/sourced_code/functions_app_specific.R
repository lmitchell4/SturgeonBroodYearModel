# Functions herein were created elsewhere and copied to here for convenience of
# use in the MitchellModelApp. Changes made to functions within this file will
# not affect the original versions of these functions (J.DuBois 10-Aug-2015)

GetStuCardLFData <- function(card_data, legit_len = 12,
                             species = c("White", "Green"), l_mit_rows = FALSE) {
  # This function processes Card data to get (essentially) a length frequency
  # distribution for each card year, raw lengths are used rather than length
  # bins (J. DuBois 14-Jul-2015)
  
  # Args:
  #   card_data:  sturgeon repord card data from 2007 to present; should contain
  #               date of capture, species, and length fields (inches and FL);
  #               should use something like dfCardDataAll
  #   legit_len:  a reasonable length (in inches) cutoff; avoids using angler-
  #               reported lengths like '3' which could be count rather than
  #               fork or total length
  #   species:    desired species for subsetting (either white or green); can
  #               use "W" or "G"
  #   l_mit_rows: row names like L. Mitchell (YYYY.length) - for compatibility
  #               with her code
  
  # Returns:
  #   a summary dataframe of length frequency distribution for each year
  
  # declare column names that are required in the dataframe (required for
  # further processing)
  needed_col_names <- c("FL_c", "Length_in", "CaptureDate", "SturgeonType")
  
  # check for column named FL_cm (if none stop)
  if (!all(needed_col_names %in% colnames(card_data))) {
    stop("Need columns ", paste0(needed_col_names, collapse = ", "),
         " in supplied dataframe.", call. = FALSE)
  }
  
  # the function YearFun is used below - so if it doesn't exist in the
  # appropriate environment, then load source code where it exists
  if (!exists("YearFun")) {
    source(file = "sourced_code/GlobalRFunctions.r",
           local = TRUE)
  }
  
  # select species type for subsetting (doing this way so the function
  # definition contains the choices for species and that default is White)
  spec_type <- match.arg(arg = species, choices = c("White", "Green"),
                         several.ok = FALSE)
  
  # index values for subsetting dataframe
  index_sub <- with(data = card_data, expr = {
    which(!is.na(Length_in) & Length_in >= legit_len &
            SturgeonType %in% spec_type)
  })
  
  # subset card data to remove length NAs and non-realistic lengths and species
  card_data_sub <- card_data[index_sub, ]
  
  # set up variable for dataframe
  Year <- YearFun(card_data_sub$CaptureDate, output = "asChar")
  ForkLength <- card_data_sub$FL_c
  
  # create dataframe of Year and ForkLength
  df_card <- data.frame(Year = Year, ForkLength = ForkLength,
                        stringsAsFactors = FALSE)
  
  # summary of fork lengths by Year and ForkLength
  df_card_sum <- plyr::ddply(.data = df_card, .variables = .(Year, ForkLength),
                             .fun = summarise,
                             Count = length(ForkLength))
  
  # convert year to numeric for use in L.Mitchell Model
  df_card_sum$Year <- as.numeric(df_card_sum$Year)
  
  # renumber for convenience (1 thru end) or like L. Mitchell
  if (l_mit_rows) {
    rownames(df_card_sum) <- paste(df_card_sum$Year,
                                   df_card_sum$ForkLength, sep = ".")
  }
  
  # function output
  df_card_sum
  #card_data_sub # for checking
}
# end GetStuCardLFData

YearFun <- function(x, output = c("asChar", "asNum")) {
  # This function extracts the year from the date field
  
  # Args:
  #   x:      the date field from which year is to be extracted
  #   output: allows user to return year as numeric or as character
  
  # Returns:
  #   four-digit year from date field either as numeric or as character
  
  # Currently no error-checking for x argument to ensure it is date
  
  output <- match.arg(arg = output, several.ok = FALSE)
  
  res <- format(x = x, format = "%Y")
  
  # Added option for numeric (05-Aug-2014)
  switch(EXPR = output,
         asChar = res <- res, # default is as.character
         asNum  = res <- as.numeric(x = res)
  )
  
  # Function output
  res  
}
# end YearFun

AssignDecade <- function(year_var) {
  # This function accepts a vector of 4-digit years (either character or
  # numeric) and returns the appropriate decade (as a string) - J. DuBois
  # 08-Jul-2015
  
  # Args:
  #   year_var: a vector of 4-digit years (either character or numeric)
  
  # Returns:
  #   the appropriate decade
  
  # for now - leaving output as character (could include option for as.factor)
  # for now - including prefix even if all years are from the same century
  #           could include code to drop prefix if years are all from same
  #           century (just an option, for now full (e.g., 1980s) will do)
  
  # ****************************************************************************
  
  # year must be 4-digit
  if (!all(nchar(year_var) == 4)) {
    stop("All years must be four-digits (YYYY).", call. = FALSE)
  }
  
  # get 2-digit prefix for each year
  prefix <- substr(x = year_var, start = 1, stop = 2)
  
  # get decade (the 3rd of the 4th digit) for each year
  decade <- substr(x = year_var, start = 3, stop = 3)
  
  # create a decades variable from "00" to "90" by 10
  decades <- structure(.Data = as.character(seq(from = 0, to = 90, by = 10)),
                       names = as.character(0:9))
  
  # given the above, I need to change the single 0 to double 0 (or "00")
  decades["0"] <- "00"
  
  # combine prefix, decades (using lookup of decade) and an 's' for output
  out <- paste0(prefix, decades[decade], "s")
  
  # function output
  out
}
# end AssignDecade

LenCat <- function(x) {
  # Added 23-Sep-2014, J. DuBois
  # This function assigns length to a catgegory based on length range or limit;
  # function based on white sturgeon total length only
  
  # Note: in R valueA < valueB (for example) will evaluate to NA if valueA is NA
  #       in the context below, if for x <= 116, x (or the fish length) is NA,
  #       then x <= 116 will evaluate to NA not "sub"; so for the way I have written
  #       the nested ifelse statement, the final 'no' (or "unk") will never evaluate
  #       have decided to keep as is for now (J. DuBois, 26-Nov-2014)
  
  # ***************************************************************************
  res <- ifelse(test = x <= 116, yes = "sub", no = 
                  ifelse(test = x %in% c(117:168), yes = "leg", no = 
                           ifelse(test = x >= 169, yes = "ovr", no = "unk"
                           )
                  )
  )
  
  # function output
  res
}
# end LenCat

