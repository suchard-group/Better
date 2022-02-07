# Jan 2022
# to run on Hoffman2

## on cluster setup
cache_dir = ''
output_dir = ''


## connection details (without using keyring)
ConnectionDetails = DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste("shinydb.cqnqzwtn5s1q.us-east-1.rds.amazonaws.com",
                 "shinydb",
                 sep = "/"),
  user = "eumaeus_readonly",
  password = "9BA3DFEE23C176",
  pathToDriver = '~/')

connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

## check if big IPC table is available in cache and get it if not saved
getIPCs(connection, 'eumaeus', cache_dir)

## run the analyses corresponding to this array id


DatabaseConnector::disconnect(connection)