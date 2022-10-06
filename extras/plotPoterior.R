# 09/06/2022
# plot posterior distributions and probs

library(wesanderson)
library(ggpubr)

source('./extras/postProcessUtils.R')
source('./extras/tryKLdivergence.R')

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
                       textSize = 18,
                        showPlot = TRUE,
                        logScale = FALSE)
## 09/13/2022 update: overlay with prior dist.
p1b = plotGBSPosteriors(allSamps,
                  adjust = adjustFlag,
                  valueRange = c(-4,4),
                  textSize = 18,
                  showPlot = TRUE,
                  logScale = FALSE,
                  showPrior = TRUE)

p2 = plotPosteriorProbs(allSamps,
                   adjust= adjustFlag,
                   colors = wes_palette("Darjeeling2"),
                   #showPlot = FALSE,
                   xpaddings = c(0.7,1),
                   textSize = 15,
                   legendPosition = c(0.8,0.5))

## 10/6/2022: add hline at 0.95 threshold
p2 = p2 + 
  geom_hline(yintercept = 0.95, 
             color = wes_palette("Royal1")[2], 
             linetype = 2, 
             size = 1)


ggarrange(p1 + rremove('xlab'),
          p2,
          nrow = 2,
          heights = c(4,2),
          align = 'h')

ggarrange(p1b + rremove('xlab'),
          p2,
          nrow = 2,
          heights = c(4,1.5))


# ## 09/07/2022
# # add KL divergence plot too (with annotations?)
# KLs = consecutiveKL(allSamps, adjust = TRUE, K = 3)
# KLs = rbind(KLs, 
#             data.frame(period = 1, KL = NA)) %>%
#   mutate(KL_label = as.character(round(KL, 3)))
# 
# p3 = ggplot(KLs, aes(x = period, y = KL)) +
#   geom_line(size = 1, color = 'gray60') +
#   geom_label(aes(label = KL_label), size = 3) + 
#   scale_x_continuous(breaks = 
#                        seq(from = min(KLs$period), 
#                            to = max(KLs$period),
#                            by = 1),
#                      expand = expansion(add = c(0.6, 1))) +
#   scale_y_continuous(expand = expansion(add = c(0.15,0.15)))+
#   labs(x='', y = 'KL div.') +
#   theme_bw()
# 
# ggarrange(p1 + rremove('xlab'), 
#           p3,
#           p2, 
#           nrow = 3,
#           heights = c(3.5,1.5,2))




