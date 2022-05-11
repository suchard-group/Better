# 05/10/2022
# helper code used for demo plots in the Bayesian sequentia testing presentation

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(wesanderson))

withAdjust = readRDS('../localCache/HCwithAdjust.rds')
withoutAdjust = readRDS('../localCache/HCwithoutAdjust.rds')

withoutAdjust = withoutAdjust %>% 
  filter(method == 'delta1=0.95 fixed',
         stats != 'delta1')

withAdjust = withAdjust %>% 
  filter(method == 'delta1=0.95 fixed',
         stats != 'delta1')

dat = rbind(withAdjust, withoutAdjust)
dat$adjust = factor(dat$adjust,
                    levels = c('without adjustment',
                               'with adjustment'))

type2cols = c(wes_palette("Zissou1")[3:4],wes_palette("Royal1")[4])
allCols = c(wes_palette("Royal1")[2], type2cols)
period_breaks = seq(from = min(dat$period_id),
                    to = max(dat$period_id),
                    by = 2)
period_labels = as.integer(period_breaks)


p = ggplot(dat, aes(x=period_id, y=y, color=stats))+
  geom_line(size = 1.5) +
  geom_point(size=2)+
  geom_hline(yintercept = 0.05,
             color = 'gray60', 
             size = 1, linetype=2)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(breaks = period_breaks, labels = period_labels)+
  labs(x='analysis period (months)', y='', color='')+
  scale_color_manual(values = allCols) +
  facet_grid(.~adjust)+
  theme_bw(base_size = 16)

print(p)