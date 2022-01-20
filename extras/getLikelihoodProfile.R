# Jan 2022
# Function: get likelihood profile from EUMAEUS saved results
# with optional function to produce a plot of the log-likelihood
# return a named double vector as the likelihood profile

## a function to get ONE likelihood profile for one specific analysis
## on one specific database
getLikelihoodProfile <- function(connection, 
                                 schema,
                                 database_id, 
                                 exposure_id, 
                                 outcome_id,
                                 analysis_id,
                                 period_id,
                                 method = "'SCCS'", # default to SCCS for our purposes
                                 plot=FALSE){
  
  sql <- "SELECT point, value 
          FROM @schema.LIKELIHOOD_PROFILE
          WHERE database_id = @database_id
          AND method = @method
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
#                            database_id = "'IBM_MDCD'", analysis_id = 1,
#                            exposure_id = 21184, outcome_id = 438945,
#                            period_id = 3, plot=TRUE)


## a function to extract multiple likelihood profiles 
## (as a list of named double vectors)
## given database




