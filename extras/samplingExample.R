# Bayesian sampling toy example

library(EvidenceSynthesis)

population <- simulatePopulations(createSimulationSettings(nSites = 1))[[1]]
cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                          data = population,
                                          modelType = "cox")
cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
likelihoodProfile <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "grid")

# evenly gridded approximate
# a double vector, names are the grid points, values are the values

# try using the simple posterior sampler here
mcmcTraces <- approximateSimplePosterior(likelihoodProfile = likelihoodProfile,
                                         priorMean = 0, priorSd = 100)

## posterior mean and 95% CI
mean(mcmcTraces$theta1)
quantile(mcmcTraces$theta1, c(0.025, 0.975))

# # try a bit on irregular likelihood profile: it works!
# mcmc2 <- approximateSimplePosterior(likelihoodProfile = lik, 
#                                     priorMean = 0, priorSd = 10)
