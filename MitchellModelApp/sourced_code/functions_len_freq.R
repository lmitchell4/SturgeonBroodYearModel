# This file houses all functions required for length frequency distribution
# analysis; functions here can be used to add a length category column (len bin)
# to any dataframe with a fish length (fork or total) field; created 11-Mar-2015
# as part of ongoing effort to clean up my code and orgranize my code - J.
# DuBois
# ******************************************************************************

# NOTE: I use length bins left included '[' but not right ')' e.g., [51,56)
# this means that this bin includes 51, 52, 53, 54, and 55 (but not 56)

# NOTE: Using delayedAssign for some of the functions below; not totally
# convinced have it dialed in, but it appears to be working as expected:
# delaying the assigment of the function until called (12-Mar-2015)

GetLenBreaks <- function(data = NULL, col_length, custom_breaks, seq.by = NULL,
                         countsNA = "always", output = "breaks") {
  # Note: 11-Mar-2015 changed name from GetBreaks to GetLenBreaks
  
  # This function calculates lengths breaks based user-supplied fish length
  # data; when data = NULL, function looks to global env for the column of
  # lengths (or col_length); breaks are set by this function and always ensure
  # the min and max values are included in a bin. In fact the lowest value of
  # the lowest bin will always be the minimum length (11-Mar-2015, J. DuBois)
  
  # assigns vector of lengths to 'len' variable
  len <- eval(expr = substitute(expr = col_length),
              envir = data, enclos = parent.frame())
  
  # gives option to set custom length breaks (for example as needed by 
  # 'WSTALKEY.xls') or if not supplied function will use min and max of fish 
  # length data; NOTE: use c(seq(from = 21, to = 186, by = 5), Inf) for white
  # sturgeon age-length key in 'WSTALKEY.xls' (J. DuBois, 12-Mar-2015)
  len_breaks <- if (missing(custom_breaks) || is.null(custom_breaks)) {
    
    # sets interval for creating bin width 
    interval <- if (missing(seq.by) || is.null(seq.by)) {
      stop("Please supply integer for 'seq.by' parameter.", call. = FALSE)
    } else {
      as.integer(seq.by)
    }
    
    # gets range (min, max) of length data
    lengthRange <- range(len, na.rm = TRUE)
    
    # setsup sequence for length breaks
    lengthBreaks <- seq(from = lengthRange[1], to = lengthRange[2], by = interval)
    
    # ensures upper bounds will always be assigned to a bin; since I am using min
    # length from range above, minimum length will always be assigned to a bin
    final_break <- if (lengthBreaks[length(lengthBreaks)] <= lengthRange[2]) {
      lengthBreaks[length(lengthBreaks)] + interval
    } else {
      NULL # max length already covered in lengthBreaks
    }
    
    # display message saying breaks will use length range
    warning("Length breaks using min and max of length data.", call. = FALSE)
    
    # creates length breaks
    c(lengthBreaks, final_break)
  } else {
    custom_breaks # use user-supplied breaks
  }
  
  # adds labels for use in AddLenCat function below; makes for a cleaner lf
  # distribution table, especially when using counts in the switch statement
  # below (added 11-Mar-2015, J. DuBois)
  break_labels <- if (!is.null(seq.by) && seq.by == 1) {
    as.character(lengthBreaks)
  } else {
    # length breaks removing last to length break removing first minus 1, 'cause
    # bins are left open right closed
    paste(len_breaks[-length(len_breaks)], len_breaks[-1] - 1, sep = "-")
    
    # below used for debugging
    #print(list(len_breaks[-length(len_breaks)], len_breaks[-1] - 1))
  }
  
  # switch function output; allows user to change output of function; select
  # 'bins' if wanting to add length bin field to dataframe
  switch(EXPR = output,
         breaks = resOut <- len_breaks,
         bins   = resOut <- AddLenCat(col_length = len, breaks = len_breaks,
                                      labels = break_labels),
         counts = resOut <- as.data.frame(
           x = table(LenBins = AddLenCat(col_length = len, breaks = len_breaks,
                                         labels = break_labels), useNA = countsNA)
           )
         )
  
  # function output
  resOut
}
# end GetLenBreaks

# Delayed functions -------------------------------------------------------
delayedAssign(x = "AddLenCat",
  function(col_length, breaks, labels) {
    # This function assigns a fish length (either in TL or FL) to a length bin;
    # the size of the bin is determined by how the user defines the 'breaks'
    # argument; AddLenCat = add length category
    
    # Bin length based on specified breaks (changed include.lowest to FALSE,
    # 11-Feb-2015, so bins were consistently the same length as the final bin
    # would have included the max value, e.g., [201-204] instead of desired
    # [201-204)
    length_cat <- cut(x = col_length, breaks = breaks, labels = labels,
                      right = FALSE, include.lowest = FALSE, ordered_result = TRUE)
    
    # added warning to alert user if any lengths were not assigned to a length bin
    if (any(!is.na(col_length) & is.na(length_cat))) {
      warning("Some lengths not included in a bin. Suggest adjusting breaks.",
              call. = FALSE)
    }
    
    # Function output
    length_cat
  }
)
# end AddLenCat

delayedAssign(x = "PlotLenFreq",
  function(data, y_var = "count", fill = "black", xlabs_odd = TRUE, ...) {
    # This function plots the length frequency distribution of fish lengths
    # using breaks from the GetLenBreaks function and ggplot; function created 
    # 11-Mar-2015 J. DuBois
    
    # capture unevaluated ... ('cause I want to pass the col_length to ggplot)
    # (from http://adv-r.had.co.nz/Computing-on-the-language.html#capturing-dots)
    dot_list <- eval(expr = substitute(alist(...)))
    
    # user-option to select count or density; not the best way to do it but for 
    # now since I'll be the user this should work pretty well (J. DuBois 
    # 11-Mar-2015)
    y <- paste0("..", y_var, "..")
    
    # establishes length breaks
    len_breaks <- GetLenBreaks(data = data, ..., output = "breaks")
    
    # deal with potential Inf value in custom length breaks; geom_histogram
    # cannot deal with Inf; assumes Inf is final value in len_breaks and
    # replaces with max length + 1 to ensure inclusion
    if (any(is.infinite(len_breaks))) {
      new_final_break <- max(data[, deparse(dot_list$col_length), drop = TRUE],
                             na.rm = TRUE) + 1
      len_breaks[which(is.infinite(len_breaks))] <- new_final_break
      
      # issue warning to let user know what is occurring with len_breaks
      warning("'Inf' length break replaced with max length value + 1.",
              call. = FALSE)     
    } else {
      len_breaks <- len_breaks
    }
    
    #NOTE: geom_histogram will throw warning if bin widths are not equal
    
    # begin steps to create plot
    lf_plot <- ggplot(data = data, mapping = aes_q(x = dot_list$col_length))
    
    # creates histogram of lengths using breaks from len_breaks
    lf_plot <- lf_plot +
      geom_histogram(mapping = aes_string(y = y), breaks = len_breaks,
                     right = FALSE, fill = fill)
    
    # sets up labels for x-axis (odds printed evens blank) - with option to use
    # len_breaks
    x_labels <- if (xlabs_odd) {
      ifelse(test = seq_along(len_breaks) %% 2 == 1,
             yes = len_breaks, no = "")
    } else {
      len_breaks
    } 
    
    # cleans up plot for display - strictly my preference
    lf_plot <- lf_plot +
      scale_x_continuous(breaks = len_breaks, labels = x_labels,
                         expand = c(0.02, 0)) +
      scale_y_continuous(expand = c(0.02, 0)) +
      theme(panel.grid.minor.x = element_blank())
    
    # function output (the plot)
    lf_plot
  }
)
# end PlotLenFreq

# delay assignment of these functions -- this did not work
#delayedAssign(x = "PlotLenFreq", value = PlotLenFreq)
