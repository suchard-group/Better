# 09/06/2022
# plot posterior distributions and probs

library(wesanderson)
library(ggpubr)

source('./extras/postProcessUtils.R')

## setups
db = 'CCAE'
me = 'HistoricalComparator'
eid = 211981
aid = 4
pr_id = 3

## 

goodOuts = c(217770473, 213966051, 217914646, 213964615, 214534072, 214211168)
outcome = goodOuts[1]

resultspath = '~/Documents/Research/betterResults/samples-shrinkMu4/'

allSamps = pullPostSamples(database_id = db,
                           method = me,
                           analysis_id = aid,
                           exposure_id = eid,
                           prior_id = pr_id,
                           outcome_id = outcome,
                           resultsPath = resultspath)

adjustFlag = TRUE

# use the GBS plotting function for now
p1 = plotGBSPosteriors(allSamps,
                  adjust = adjustFlag,
                  valueRange = c(-4,4),
                  showPlot = FALSE,
                  logScale = FALSE)

p2 = plotPosteriorProbs(allSamps,
                   adjust= adjustFlag,
                   colors = wes_palette("Darjeeling2"),
                   #showPlot = FALSE,
                   xpaddings = c(0.5,1))

ggarrange(p1 + rremove('xlab'), 
          p2, 
          nrow = 2,
          heights = c(3,1.5))
