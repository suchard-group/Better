# Jan 2022
# Function: get likelihood profile from EUMAEUS saved results
# with optional function to produce a plot of the log-likelihood
# return a named double vector as the likelihood profile

library(tidyverse)
library(stringr)

## a function to get ONE likelihood profile for one specific analysis
## on one specific database
## !!returns numeric(0) if no match exists in the likelihood profile table
getLikelihoodProfile <- function(connection, 
                                 schema,
                                 database_id, 
                                 exposure_id, 
                                 outcome_id,
                                 analysis_id,
                                 period_id,
                                 method = "SCCS", # default to SCCS for our purposes
                                 plot=FALSE){
  
  sql <- "SELECT point, value 
          FROM @schema.LIKELIHOOD_PROFILE
          WHERE database_id = '@database_id'
          AND method = '@method'
          AND analysis_id = @analysis_id
          AND exposure_id = @exposure_id
          AND outcome_id = @outcome_id
          AND period_id = @period_id"
  
  ## using the updated "render" & "translate" functions instead
  ## they return chatacter strings directly; no need to extract from a list
  sql <- SqlRender::render(sql, 
                           schema = schema,
                           database_id = database_id,
                           method = method,
                           analysis_id = analysis_id,
                           exposure_id = exposure_id,
                           outcome_id = outcome_id,
                           period_id = period_id)
  sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
  lik <- DatabaseConnector::querySql(connection, sql)
  
  # check to see if result is empty
  if(nrow(lik) == 0){
    res = numeric(0)
  }else{
    # split string and convert to numbers
    points = str_split(lik$POINT,';') %>% 
      unlist() %>% as.numeric()
    values = str_split(lik$VALUE,';') %>% 
      unlist() %>% as.numeric()
    
    # # create a data frame as result
    # res = data.frame(point = points, value = values)
    
    # UPDATED (for BEAST usage)
    # create a double vector as result
    res = values
    names(res) = as.character(points)
    
    # plot if...
    if(plot){
      # construct dataframe for plotting
      res.dat = data.frame(point = points, value = values)
      g = ggplot(res.dat, aes(x=point, y=value)) +
        geom_line(size = 0.8) +
        labs(x = 'parameter value', y = 'log-likelihood') +
        theme_bw(base_size=14)
      print(g)
    }
  }
  
  # return the likelihood profile as a double vector
  res
}

# example
# lik = getLikelihoodProfile(connection, "eumaeus",
#                            database_id = "IBM_MDCD", analysis_id = 1,
#                            exposure_id = 21184, outcome_id = 438945,
#                            period_id = 3, plot=TRUE)


## a function to extract multiple likelihood profiles 
## (as a list of named double vectors)
## given database_id, method, analysis, period and exposure
## optional: can subset on outcomes
## process: if to split the string points and values and return a list
##          (FALSE: return the raw queried table w/o splitting the strings)
## Feb 2022 update: can allow (1 period, multi analyses) OR (1 analysis, multi periods)

# April 2023 update: add a cleanUp functionality to filter out NAs in `POINT` column
getMultiLikelihoodProfiles <- function(connection, 
                                       schema,
                                       database_id,
                                       exposure_id, 
                                       analysis_id = NULL,
                                       period_id = NULL,
                                       outcome_ids = NULL,
                                       method = "SCCS",
                                       process = FALSE,
                                       cleanUp = TRUE){
  # query all likelihood profiles needed
    ## if multiple analyses...
  if(is.null(analysis_id)){
    sql <- "SELECT * 
          FROM @schema.LIKELIHOOD_PROFILE
          WHERE database_id = '@database_id'
          AND method = '@method'
          AND exposure_id = @exposure_id
          AND period_id = @period_id"
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = database_id,
                             method = method,
                             exposure_id = exposure_id,
                             period_id = period_id)
  }else{
    ## otherwise, multiple periods
    sql <- "SELECT * 
          FROM @schema.LIKELIHOOD_PROFILE
          WHERE database_id = '@database_id'
          AND method = '@method'
          AND exposure_id = @exposure_id
          AND analysis_id = @analysis_id"
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = database_id,
                             method = method,
                             exposure_id = exposure_id,
                             analysis_id = analysis_id)
  }
  
  LPs = DatabaseConnector::querySql(connection, sql)
  cat('Likelihood profiles extracted.\n')
  
  # filter on needed subsets...
  
  # if(!is.null(exposure_ids)){
  #   LPs  = LPs %>% filter(exposure_id %in% exposure_ids)
  # }
  if(!is.null(outcome_ids)){
    LPs  = LPs %>% filter(OUTCOME_ID %in% outcome_ids)
  }
  # if(!is.null(analysis_ids)){
  #   LPs  = LPs %>% filter(analysis_id %in% analysis_ids)
  # }
  # if(!is.null(period_ids)){
  #   LPs  = LPs %>% filter(period_id %in% period_ids)
  # }
  
  # clean it up and organize to a list
  # if(process){
  #   # don't do anything for now
  #   res = LPs
  # }else{
  #   res = LPs
  # }
  
  # clean it up; filter out "NA" entries in the POINT column
  if(cleanUp){
    if(nrow(LPs) == 0){
      return(LPs)
    }
    
    bad_indices = grep("*NA*", LPs$POINT)
    if(length(bad_indices) > 0){
      LPs = LPs[-bad_indices,]
    }
  }
  
  # lower case the col names
  #names(res) = tolower(names(res))
  LPs
  
}

## try it and time it
## it seems more reasonable now
# t1 = Sys.time()
# LPs = getMultiLikelihoodProfiles(connection, 'eumaeus',
#                                  database_id = "IBM_MDCD",
#                                  exposure_id = 21184, # H1N1 vaccine
#                                  method = "SCCS",
#                                  period_id = NULL,
#                                  analysis_id = 1)
# Sys.time() - t1

# July 2022 update: post-process multiple likelihood profiles 
# and make them into a list of dataframes (by outcome_ids)
postProcessLPs <- function(LPs, 
                           period_ids = NULL,
                           analysis_ids = NULL,
                           name_by_outcome = TRUE){
  
  names(LPs) = tolower(names(LPs))
  if(!is.null(analysis_ids)){
    LPs = LPs %>% filter(analysis_id %in% analysis_ids)
  }
  if(!is.null(period_ids)){
    LPs = LPs %>% filter(period_id %in% period_ids)
  }
  
  LPs = LPs %>% arrange(outcome_id)
  LPs$index = c(1:nrow(LPs))
  
  LPlist = split(LPs, LPs$index)
  for(i in 1:length(LPlist)){
    points = str_split(LPlist[[i]]$point,';') %>% 
      unlist() %>% as.numeric()
    values = str_split(LPlist[[i]]$value,';') %>% 
      unlist() %>% as.numeric()
    LPlist[[i]] = data.frame(point = points, value = values)
  }
  
  if(name_by_outcome){
    names(LPlist) = as.character(LPs$outcome_id)
  }
  
  LPlist
}


## function to get a particular likelihood profile from the pulled raw datatable
selectLikelihoodProfileEntry <- function(df,
                                         database_id,
                                         method,
                                         exposure_id, 
                                         period_id,
                                         outcome_id,
                                         analysis_id,
                                         plot = FALSE){
  
  # check if df has any rows at all
  if(nrow(df) == 0){
    # returns empty double vector if nothing
    res = numeric(0)
  }else{
    lik = df %>% 
      filter(DATABASE_ID == database_id, 
             METHOD == method,
             EXPOSURE_ID == exposure_id,
             PERIOD_ID == period_id,
             OUTCOME_ID == outcome_id, 
             ANALYSIS_ID == analysis_id)
    if(nrow(df) == 0){
      # return empty if no results
      res = numeric(0)
    }else{
      # split string and convert to numbers
      ## report rows we've got here
      #cat('rows returned: ', nrow(lik),'\n')
      points = str_split(lik$POINT,';') %>% 
        unlist() %>% as.numeric()
      values = str_split(lik$VALUE,';') %>% 
        unlist() %>% as.numeric()

      # create a double vector as result
      res = values
      names(res) = as.character(points)
      
      if(plot){
        # construct dataframe for plotting
        res.dat = data.frame(point = points, value = values)
        g = ggplot(res.dat, aes(x=point, y=value)) +
          geom_line(size = 0.8) +
          labs(x = 'parameter value', y = 'log-likelihood') +
          theme_bw(base_size=14)
        print(g)
      }
      
    }
  }
  # return result (either has empty or named double vector)
  res
}

# # try it
# lik = selectLikelihoodProfileEntry(LPs, 'IBM_MDCD', 'SCCS',
#                                    21184, 9,
#                                    outcome_id = 10109,
#                                    analysis_id = 1, plot=TRUE)
