# 11/21/2023
# plot exposure counts and aggregated exposure time

library(dplyr)
library(ggplot2)

# some test code first----

## setup
fpath = '~/Documents/Research/better_gbs/Results_CCAE/estimate.csv'

database_id = 'CCAE'
method = 'HistoricalComparator'
analysis_id = 4
exposure_id = 211983
outcome_id = 343

cachePath = './localCache/'

exposure_select = c(21184, 21185, 211983, 211833) 
### these have reletively same overall counts...


## parse dataframe
estimates = read_csv(fpath)
names(estimates) = SqlRender::camelCaseToSnakeCase(names(estimates))

this.est = estimates %>% 
  filter(database_id == !!database_id,
         method == !!method,
         analysis_id == !!analysis_id,
         exposure_id %in% exposure_select,
         outcome_id == !!outcome_id) %>%
  select(period_id, exposure_id, 
         exposure_subjects, exposure_days)

exposures = readRDS(file.path(cachePath, 'exposures.rds'))

this.est = this.est %>% 
  left_join(exposures, by = 'exposure_id') %>%
  select(period_id, exposure_id, 
         exposure_subjects, exposure_days,
         exposure = base_exposure_name)

## plot
xbreaks = seq(from = 0, to = max(this.est$period_id), by = 3)
theColors = c(wes_palette("Darjeeling1")[2:3], 
              wes_palette("Darjeeling2")[2],
              wes_palette("GrandBudapest1")[2])

count_labels = this.est %>% 
  group_by(exposure) %>%
  summarise(period = max(period_id) - 1, 
            counts = max(exposure_subjects) + 2*1e4) %>%
  ungroup()

ggplot(this.est, aes(x = period_id, 
                     y = exposure_subjects, 
                     color = exposure)) +
  geom_line(size = 1.2) +
  # geom_label(data = count_labels, color = 'black', 
  #            size = 4,
  #            aes(x=period, y = counts, label = exposure)) +
  scale_x_continuous(breaks = xbreaks, limits = c(1, max(this.est$period_id))) +
  guides(color=guide_legend(nrow=2,byrow=FALSE))+
  labs(x='Analysis period', y = 'Exposure subject count',
       color = 'Vaccine exposure') +
  theme_bw(base_size = 14) +
  scale_color_manual(values = theColors) +
  theme(legend.position = 'bottom')

ggplot(this.est, aes(x = period_id, 
                     y = exposure_days, 
                     color = exposure)) +
  geom_line(size = 1.2) +
  # geom_label(data = count_labels, color = 'black', 
  #            size = 4,
  #            aes(x=period, y = counts, label = exposure)) +
  scale_x_continuous(breaks = xbreaks, limits = c(1, max(this.est$period_id))) +
  guides(color=guide_legend(nrow=2,byrow=FALSE))+
  labs(x='Analysis period', y = 'Aggregated exposure days',
       color = 'Vaccine exposure') +
  theme_bw(base_size = 14) +
  scale_color_manual(values = theColors) +
  theme(legend.position = 'bottom')
  