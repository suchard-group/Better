# 11/21/2023
# plot exposure counts and aggregated exposure time

library(dplyr)
library(ggplot2)
library(readr)

# some test code first----

## setup
fpath = '~/Documents/Research/better_gbs/Results_CCAE/estimate.csv'
fpath = '~/Documents/Research/better_gbs/Results_MDCD/estimate.csv'

#database_id = 'CCAE'
method = 'HistoricalComparator'
analysis_id = 1
#exposure_id = 211983
outcome_id = 343

cachePath = './localCache/'

exposure_select = c(21184, 21185, 21214, 21215, 211983, 211833) 
theColors = c(wes_palette("Darjeeling1")[2:3], 
              wes_palette("Darjeeling2")[2],
              wes_palette("GrandBudapest1")[2],
              wes_palette("IsleofDogs1")[1:2])
### these have reletively same overall counts...


## parse dataframe
estimates = read_csv(fpath)
names(estimates) = SqlRender::camelCaseToSnakeCase(names(estimates))

this.est = estimates %>% 
  filter(method == !!method,
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
  scale_y_log10(limits = c(1, 5*1e6))+
  guides(color=guide_legend(nrow=3,byrow=TRUE))+
  labs(x='Analysis period (in months)', y = 'Exposure subject count',
       color = 'Vaccine exposure', title = sprintf('Database: %s',database_id)) +
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


## 02/02/2023----
# plot all databases' subject counts over time in one big plot
library(gridExtra)

plist = list()
database_list = c('CCAE', 'MDCR', 'MDCD', 'OptumEHR', 'OptumDod')

overall_legend = NULL

exposures = readRDS(file.path(cachePath, 'exposures.rds'))

for(db in database_list){
  fpath = sprintf("~/Documents/Research/better_gbs/Results_%s/estimate.csv", db)
  
  estimates = read_csv(fpath)
  names(estimates) = SqlRender::camelCaseToSnakeCase(names(estimates))
  
  this.est = estimates %>% 
    filter(method == !!method,
           analysis_id == !!analysis_id,
           exposure_id %in% exposure_select,
           outcome_id == !!outcome_id) %>%
    select(period_id, exposure_id, 
           exposure_subjects, exposure_days)
  
  this.est = this.est %>% 
    left_join(exposures, by = 'exposure_id') %>%
    select(period_id, exposure_id, 
           exposure_subjects, exposure_days,
           exposure = base_exposure_name)
  
  ## plot
  xbreaks = seq(from = 0, to = max(this.est$period_id), by = 3)
  
  this.plot = 
    ggplot(this.est, aes(x = period_id, 
                       y = exposure_subjects, 
                       color = exposure)) +
    geom_line(size = 1.2) +
    # geom_label(data = count_labels, color = 'black', 
    #            size = 4,
    #            aes(x=period, y = counts, label = exposure)) +
    scale_x_continuous(breaks = xbreaks, limits = c(1, max(this.est$period_id))) +
    scale_y_log10(limits = c(1, 5*1e6))+
    guides(color=guide_legend(nrow=6,byrow=TRUE))+
    labs(x='Analysis period (in months)', y = 'Number of\nexposure subjects',
         color = 'Vaccine exposure', title = sprintf('Database: %s',db)) +
    theme_bw(base_size = 14) +
    scale_color_manual(values = theColors)
  
  if(which(database_list == db) == 1){
    overall_legend = lemon::g_legend(this.plot+
                                       theme(legend.position = 'right',
                                             legend.text=element_text(size=12)))
  }

  this.plot = this.plot + theme(legend.position = 'hidden')
  
  plist[[which(database_list == db)]] = this.plot
  
}

plist[[length(plist)+1]] = overall_legend
grid.arrange(grobs = plist,
             ncol = 2)


# 04/13/2023: use data characterization results local copy & add CUIMC -----
# still, one bit plot for exposure counts over time

estimates = readRDS('./localCache/negativeControlEstimatesCharacteristics-2.rds')

library(gridExtra)

exposures = readRDS(file.path('./localCache/', 'exposures.rds'))
theColors = c(wes_palette("Darjeeling1")[2:3], 
              wes_palette("Darjeeling2")[2],
              wes_palette("GrandBudapest1")[2],
              wes_palette("IsleofDogs1")[1:2])

exposed_counts = estimates %>% 
  left_join(exposures, by = 'exposure_id') %>%
  select(database_id, period_id, exposure_id, 
         exposure_subjects, exposure_days,
         exposure = base_exposure_name) %>% 
  group_by(database_id, period_id, exposure) %>%
  summarise(exposure_count = max(exposure_subjects)) %>%
  ungroup()

plist = list()
overall_legend = NULL

all_databases = unique(exposed_counts$database_id)
for(db in all_databases){
  this.exposed = exposed_counts %>% filter(database_id == db)
  
  if(stringr::str_starts(db,'IBM')){
    dbname = unlist(stringr::str_split(db, '_'))[2]
  }else{
    dbname = db
  }
  
  # 04/24/2023: change "Optum" name to "Clinformatics" instead
  if (dbname == 'OptumDod'){
    dbname = "Clinformatics"
  }
  
  this.p = ggplot(this.exposed, aes(x = period_id, 
                                y = exposure_count, 
                                color = exposure)) +
    geom_line(size = 1.2) +
    scale_x_continuous(breaks = seq(from=3, to=12, by=3), limits = c(1, 12)) +
    scale_y_log10(limits = c(1, 5*1e6))+
    guides(color=guide_legend(nrow=6,byrow=TRUE))+
    labs(x='Analysis period (in months)', y = 'Number of\nexposed subjects',
         color = 'Vaccine exposure', title = sprintf('Database: %s',dbname)) +
    theme_bw(base_size = 16) +
    scale_color_manual(values = theColors)
  
  if(db == 'CCAE'){
    overall_legend = lemon::g_legend(this.p+
                                       theme(legend.position = 'right',
                                             legend.text=element_text(size=12))+
                                       guides(color=guide_legend(nrow=2,byrow=TRUE)))
  }
  this.p = this.p + theme(legend.position = 'hidden')
  plist[[which(all_databases == db)]] = this.p
}


p_main = grid.arrange(grobs = plist,
                      ncol = 2)

p_all = grid.arrange(p_main, overall_legend, 
             nrow = 2,
             heights = c(8,1))

ggsave("~/Documents/Research/betterResults/plots/exposureCounts.pdf",
       p_all,
       height = 10, width = 10)

  