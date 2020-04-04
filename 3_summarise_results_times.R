# 3_summarise_results_times.R
# summarise the results of the microsimulations
# examining transition times
# March 2020
library(ggplot2)
library(stringr)
library(dplyr)

## create times between transitions
# first make data frame of starting (initial) times
initial = expand.grid(ID=1:meta$N.start, simnum=1:max(transitions$simnum)) # all combinations
initial = mutate(initial, day=-1, transitionTime=as.numeric(as.Date(meta$first.date, '%d/%m/%Y')))
# calculate times
time.between = bind_rows(initial, transitions) %>%
  group_by(simnum, ID) %>%
  arrange(simnum, ID, day, transitionTime) %>%
  mutate(journey = paste(From, '->', To, sep=''),
         lag = lag(transitionTime),
         diff = transitionTime - lag + 0.5) %>% # add half a day - transitions on same day
  filter(day > -1) %>% # remove starting states as no difference to calculate
  ungroup()


## to fix
# histograms of times - arrange in single column
hplot = ggplot(filter(time.between, From=='I3'), aes(x=diff))+
  geom_histogram(bins=20, fill='sky blue')+
  xlab('Time to transition (days)')+
  ylab('Count')+
  facet_wrap(~journey, scales='free_y')+
  theme_bw()
hplot
jpeg('figures/FromI3.jpg', width=5, height=4, units='in', res=300)
print(hplot)
dev.off()


## summary stats on transition times
# average per simulation
sim.times = group_by(time.between, simnum, From, To) %>%
  summarise(n=n(), meant=mean(diff), mediant = median(diff), Q1t=quantile(diff, 0.25), Q3t=quantile(diff, 0.75))
group_by(sim.times, From, To) %>%
  summarise(n=n(), mean=mean(meant))
# plot means
mplot = ggplot(data=filter(sim.times, !From %in% c('HR','HS','HE','H1')), aes(x=factor(From), y=meant, fill=factor(To)))+
  geom_boxplot()+
  ylab('Mean time in days')+
  xlab('From')
mplot


#
hplot = ggplot(filter(sim.times, From=='I3'), aes(meant))+
  geom_histogram()+
  facet_wrap(~To)
hplot


### randomly plot one simulation ###
# histograms of times - arrange in single column
random = time.between[time.between$simnum == 6,]
hplot = ggplot(filter(random, From=='I2'), aes(x=diff))+
  geom_histogram(bins=20, fill='sky blue')+
  xlab('Time to transition (days)')+
  ylab('Count')+
  facet_wrap(~To, scales='free_y')+
  theme_bw()
hplot
