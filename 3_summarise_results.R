# 3_summarise_results.R
# summarise the results of the microsimulation
# March 2020
library(ggplot2)
library(stringr)
library(dplyr)
# get the data (from 0_microsim_covid)
load(file='data/all_results.RData') # from 2_combine_simulations.R
load(file='data/all_results_vary.RData') # from 2_combine_simulations.R, verion with plus/minus 10% in parameters
# ODE results from Alison's model, from runSpread_compare.R
load('data/AlisonResults.RData')

# remove genders to examine overall transitions
transitions = mutate(all_results, 
                     From = str_remove(From, pattern='f/|m/'), # remove genders
                     To = str_remove(To, pattern='f/|m/')) 

# counts per Simulation of transitions
counts = group_by(transitions, simnum, From, To) %>%
  summarise(n=n()) %>%
  arrange(simnum, n)
# summary stats on counts across simulations (inter-quartile range and 95% CI)
count.summary = group_by(counts, From, To) %>%
  summarise(median=median(n), Q1=quantile(n, 0.25), Q3=quantile(n, 0.75), lower=quantile(n, 0.05), upper=quantile(n, 0.95)) %>%
  arrange(-median) 
#%>%   filter(To != 'dead') # do not worry about long-term deaths
count.summary

# counts per Simulation of numbers travelling through each state
count.final = group_by(transitions, simnum, To) %>%
  summarise(n=n()) %>%
  arrange(simnum, n)
# summary stats on counts across simulations (inter-quartile range and 95% CI)
count.summary.final = group_by(count.final, To) %>%
  summarise(median=median(n), Q1=quantile(n, 0.25), Q3=quantile(n, 0.75), lower=quantile(n, 0.05), upper=quantile(n, 0.95)) %>%
  arrange(-median)
count.summary.final

## table of transitions across all simulations
with(transitions, table(From, To))

## times between transitions
# first make data frame of starting (initial) times
initial = expand.grid(ID=1:meta$N.start, simnum=1:max(transitions$simnum)) # all combinations
initial = mutate(initial, day=-1, transitionTime=as.numeric(as.Date(meta$first.date, '%d/%m/%Y')))
# calculate times
time.between = bind_rows(initial, transitions) %>%
  group_by(ID) %>%
  arrange(ID, day, transitionTime) %>%
  mutate(lag = lag(transitionTime),
         diff = transitionTime - lag) %>%
  filter(day > -1) %>% # remove starting states as no difference to calculate
  ungroup()

# calculate median times per simulation and plot that 

## to fix
# histograms of times - arrange in single column
hplot = ggplot(filter(time.between, To=='R'), aes(x=diff))+
  geom_histogram(bins=20)+
  xlab('Time to transition (days)')+
  facet_wrap(~From, scales='free_y')+
  theme_bw()
hplot
# summary stats on time
group_by(time.between, From, To) %>%
  summarise(n=n(), median = median(diff), Q1=quantile(diff, 0.25), Q3=quantile(diff, 0.75))

## plot transitions over time
# convert to cumulative across days
times = group_by(transitions, simnum, To, transitionTime) %>%
  summarise(n=n()) %>%
  arrange(simnum, To, transitionTime) %>%
  group_by(simnum, To) %>%
  mutate(cum = cumsum(n),
         day = as.numeric(transitionTime - min(transitionTime)),
         ToNice = factor(To, levels=c('E','I1','I2','I3','R','D','dead'), # make nicer labels
                         labels = c('Exposed','Mild infection','Severe infection','Critical infection','Recovered','Died - COVID','Died - Other'))) %>%
  ungroup()
# a) cumulative numbers
colours = grey(runif(50, 0.2, 0.8))
cplot = ggplot(data=times, aes(x=day, y=cum, group=factor(simnum), col=factor(simnum)))+
  geom_step(size=0.5)+
  scale_color_manual(NULL, values=colours)+
  xlab('Day')+
  ylab('Cumulative number')+
  coord_cartesian(xlim=c(0,150))+
  facet_wrap(~ToNice, scales='free_y')+
  theme_bw()+
  theme(legend.position = 'none')
cplot
#
jpeg('figures/microsim.numbers.jpg', width=5, height=4, units='in', res=300)
print(cplot)
dev.off()

# b) daily incidence numbers (not volume)
dplot = ggplot(data=times, aes(x=day, y=n, group=factor(simnum)))+
  geom_step(col='grey', size=0.5)+
  xlab('Day')+
  ylab('Number')+
  facet_wrap(~ToNice, scales='free_y')+
  theme_bw()+
  theme(legend.position = 'none')
dplot

# age at death for other deaths
group_by(transitions, To) %>%
  summarise(mean = mean(transitionAge), Q1=quantile(transitionAge, 0.25), Q3=quantile(transitionAge, 0.75))

# TO DO, recreate timelines per person to get density of time in each state.


# Plot with Alison's numbers on top
recovered.microsim = filter(times, To=='R') %>% 
  mutate(source='Microsimulation') %>%
  select(-transitionTime, -To, -ToNice) %>%
  dplyr::rename('time' = 'day',
         'value'='cum') # to match Alison's labels
  
recovered.ode = filter(out, variable=='R') %>%
  mutate(source='ODE',
         simnum = 51 # for line colouring
         )
to.plot = bind_rows(recovered.microsim, recovered.ode)
colours = c(grey(runif(50, 0.2, 0.8)),'dark blue')
compare.plot = ggplot(data=to.plot, aes(x=time, y=value, group=factor(simnum), col=factor(simnum)))+
  geom_step(size=0.5)+
  xlab('Day')+
  scale_color_manual(NULL, values=colours)+
  ylab('Cumulative number')+
  coord_cartesian(xlim=c(0,150))+
  theme_bw()+
  theme(legend.position = 'none')
compare.plot
