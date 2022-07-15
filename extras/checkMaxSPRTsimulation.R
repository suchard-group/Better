# 07/14/2022
# simulation to check error rates of MaxSPRT

library(Sequential)

# use a Poisson model first

## a function to simulate Poisson data and run sequential test
simulatePoissonTest <- function(S, numLooks=12, 
                                Nbatch=100, 
                                expectedCount = 5,
                                effectSize = 1, 
                                alphaSpend = 'Wald',
                                cachePath = './localCache/'){
  
  rejects = 0
  
  for(s in 1:S){
    cat(sprintf('Simulation round %s:...\n', s))
    testName = sprintf('Poisson%s', s)
    
    if(file.exists(file.path(cachepath, sprintf('%s.txt',testName)))){
      file.remove(file.path(cachepath, sprintf('%s.txt',testName)))
    }
    
    AnalyzeSetUp.Poisson(name = testName, 
                         SampleSize = expectedCount * numLooks,
                         M = 4, # recommended by Kulldorff and Silva
                         AlphaSpendType = alphaSpend,
                         address = cachePath)
    
    #AllRes = NULL
    
    for(i in 1:numLooks){
      ## simulate Poisson count data
      obsCounts = sum(rpois(Nbatch, expectedCount/Nbatch))
      
      ## test
      res = Analyze.Poisson('Poisson1', test = i, 
                            mu0 = expectedCount, 
                            events = obsCounts)
      
      if(res[1,'Reject H0'] == 'Yes'){
        cat(sprintf('Null rejected at %sth data look.\n', i))
        rejects = rejects + 1
        break
      }
    }
    
    rm()
  }
  
  return(rejects)
}


## try simulation run

cachepath = '~/Documents/Research/better/localCache/'

rejects = simulatePoissonTest(1, numLooks = 12, 
                              Nbatch = 100,
                              expectedCount = 5,
                              effectSize = 1, 
                              alphaSpend = 'Wald',
                              cachePath = cachepath)
# somehow there is an internal stuff that doesn't allow me to repeat testing inside the function???

S = 100
allRejects = 0
for(s in 1:S){
  reject = simulatePoissonTest(1, numLooks = 12, 
                               Nbatch = 100,
                               expectedCount = 5,
                               effectSize = 1, 
                               alphaSpend = 'Wald',
                               cachePath = cachepath)
  allRejects = allRejects + reject
}

Type1error = allRejects/S
