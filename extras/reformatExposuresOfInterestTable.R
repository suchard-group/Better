# transform date format for exposures-of-interest table...

eois  = readr::read_csv(system.file("settings", "ExposuresOfInterest.csv", package = "Better"), 
                        col_types = readr::cols()) %>% 
  select(exposureName, startDate, endDate, historyStartDate, historyEndDate)%>%
  mutate(startDate = format(strptime(startDate, format = '%d-%m-%Y'), "%m/%d/%Y"),
         endDate = format(strptime(endDate, format = '%d-%m-%Y'), "%m/%d/%Y"),
         historyStartDate = format(strptime(historyStartDate, format = '%d-%m-%Y'), "%m/%d/%Y"),
         historyEndDate = format(strptime(historyEndDate, format = '%d-%m-%Y'), "%m/%d/%Y"))

colnames(eois) <- SqlRender::camelCaseToTitleCase(colnames(eois))


print(xtable::xtable(eois), include.rownames=FALSE)
