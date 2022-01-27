# Jan 2022
# a test script to run stuff on Hoffman2

library(EvidenceSynthesis)
source('./extras/getLikelihoodProfile.R')
source('./extras/fitNegativeControlDistribution.R')

library(foreach)
library(doParallel)
registerDoParallel()

library(ParallelLogger)

## 0. get array_id
array_id = Sys.getenv('SGE_TASK_ID')
cat('Array ID is', array_id, '\n\n')

## 1. try connecting to EUMAEUS results server
ConnectionDetails = DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste("shinydb.cqnqzwtn5s1q.us-east-1.rds.amazonaws.com",
                 "shinydb",
                 sep = "/"),
  user = "eumaeus_readonly",
  password = "9BA3DFEE23C176",
  pathToDriver = '~/')

connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

sql = "SELECT COUNT(*) from eumaeus.ANALYSIS"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
num.analysis <- DatabaseConnector::querySql(connection, sql)
cat('Total number of EUMAEUS analyses: ', num.analysis$COUNT, '\n\n')

DatabaseConnector::disconnect(connection)

## 2. try run foreach
xx = foreach(i=1:3, .combine = 'c') %dopar% {
  rnorm(5)
}

cat(sprintf('Total length of combined vector is %s, with mean %.3f\n\n', 
            length(xx), mean(xx)))


## 3. try using ParallelLogger
fun <- function(x) {
  x^2
}
cluster <- makeCluster(numberOfThreads = 2)
yy = unlist(
  clusterApply(cluster, 1:10, fun)
)
cat('ParallelLogger test results: y=', paste(yy, collapse = ','),'.\n')
stopCluster(cluster)
