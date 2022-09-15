# some illustrative plots for BETTER presentation at OHDSI Symposium

library(wesanderson)

# 1. plot for planned vs actual Type 1 error rate

dat = data.frame(x=as.factor(1:2), y = c(0.05, 0.25))

ggplot(dat, aes(x=x, y=y, fill = x)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0.05, color='gray50', linetype=2, size = 1)+
  scale_y_continuous(limits = c(0,0.35))+
  scale_x_discrete(labels = NULL, breaks = NULL)+
  scale_fill_manual(values = wes_palette("Royal1")[c(2,4)])+
  labs(x='',y='')+
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')

# 2. bias adjustment diagram

fakeNCs = rnorm(90, mean = 1.2, sd = 1)
dat = data.frame(nc = fakeNCs)
ggplot(dat, aes(x=nc))+
  geom_density(fill = wes_palette("Royal2")[3], alpha = 0.3) +
  geom_point(aes(y = 0.03), shape = 5)+
  scale_x_continuous(limits = c(-2.5,2.5))+
  scale_y_continuous(breaks = NULL)+
  labs(x='',y='')+
  theme_minimal(base_size = 15)

fakebetas = rnorm(5000, mean = 2.7, sd = 1.5)
dat2 = data.frame(beta = fakebetas)
ggplot(dat2, aes(x=beta))+
  geom_density(fill = wes_palette("Royal1")[4], alpha = 0.3) +
  #geom_point(aes(y = 0.03), shape = 5)+
  scale_x_continuous(limits = c(-2.5,5))+
  scale_y_continuous(breaks = NULL)+
  labs(x='',y='')+
  theme_minimal(base_size = 15)

adjbetas = fakebetas - rnorm(5000, mean = 1.2, sd = 1)
dat3 = data.frame(beta = adjbetas)
ggplot(dat3, aes(x=beta))+
  geom_density(fill = wes_palette("Royal1")[2], alpha = 0.3) +
  #geom_point(aes(y = 0.03), shape = 5)+
  scale_x_continuous(limits = c(-2.5,5))+
  scale_y_continuous(breaks = NULL)+
  labs(x='',y='')+
  theme_minimal(base_size = 15)

