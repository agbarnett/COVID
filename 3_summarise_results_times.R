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
initial = mutate(initial, 
                 day = -1, 
                 transitionTime = as.numeric(as.Date(meta$first.date, '%d/%m/%Y')), 
                 From = ifelse(ID <= meta$starting.numbers, 'E', 'S')) # first numbers are exposed, others susceptible

# calculate times
time.between = bind_rows(initial, transitions) %>%
  group_by(simnum, ID) %>%
  arrange(simnum, ID, day, transitionTime) %>%
  mutate(journey = paste(From, '->', To, sep=''),
         lag = lag(transitionTime),
         diff = transitionTime - lag + 0.5) %>% # add half a day - transitions on same day
  filter(day > -1) %>% # remove starting states as no difference to calculate
  ungroup()

### Section 1: times between states ###

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

### Section 2: incident numbers in each state over time ###

## create row for each person on each day and their state
# expand all observed days for each person
all.days = bind_rows(initial, transitions) %>%
  group_by(simnum, ID) %>%
  expand(days = full_seq(day, 1)) %>% # takes a while
  ungroup()
# make a smaller data set of just states moved to
smaller = bind_rows(initial, transitions) %>%
  mutate(To = ifelse(is.na(To)==TRUE, From, To),
         day = ifelse(day==-1, 0, day)) %>% # move initial transition to day zero
  select('simnum','ID','day','To') %>%
  rename('days' = 'day')
# expand over days and count people in each state on each day
to.plot = left_join(all.days, smaller, by=c('simnum','ID','days')) %>%
  fill(simnum, ID, days, To) %>% #fill in missing states
  filter(days >=0) %>%
  group_by(simnum, days, To) %>%
  summarise(count = n()) %>%
  rename('State' = 'To')
# plot numbers over time  
colours = grey(runif(50, 0.2, 0.8))
to.plot = filter(to.plot, !State %in%c('H1','HS','HR','HE'))
dplot = ggplot(data=to.plot, aes(x=days, y=count, col=factor(simnum)))+
  geom_step()+
  scale_color_manual(NULL, values=colours)+
  facet_wrap(~State, scales = 'free_y')+
  theme(legend.position = 'none')
dplot
