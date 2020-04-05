# 3_summarise_results.R
# summarise the results of the microsimulations
# March 2020
library(ggplot2)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999") # colours
library(stringr)
library(dplyr)
library(tidyr)
# get the data (from 1_microsim_covid[...] on lyra)
load(file='data/all_results.RData') # from 2_combine_simulations.R
load(file='data/all_results_vary.RData') # from 2_combine_simulations.R, verion with plus/minus 10% in parameters

# remove genders to examine overall transitions
transitions = mutate(all_results, 
                     From = str_remove(From, pattern='f/|m/'), # remove genders
                     To = str_remove(To, pattern='f/|m/')) 

### Section 0: seasonal pattern ###

days = as.Date(meta$first.date,  '%d/%m/%Y') + 0:meta$max.day
seasonal = 1 + meta$seas.amp*cos(2*pi*(as.numeric(days) -  as.numeric(as.Date(meta$seas.phase))) /365.25)
plot(days, seasonal, type='o')

### Section 1: counts ###

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

### Section 2: plot transitions over time ###

# convert to cumulative across days
times = group_by(transitions, simnum, To, transitionTime) %>%
  filter(To !='S') %>% # not interested in returns to susceptible
  summarise(n=n()) %>%
  arrange(simnum, To, transitionTime) %>%
  group_by(simnum, To) %>%
  mutate(cum = cumsum(n),
         day = as.numeric(transitionTime - min(transitionTime)),
         ToNice = factor(To, levels=c('E','I1','I2','I3','R','D','HR','HS','H1','HE','dead'), # make nicer labels
                         labels = c('Exposed','Mild infection','Severe infection','Critical infection','Recovered','Died - COVID','Hospital R','Hospital S','Hospital I1','Hospital E','Died - Other'))) %>%
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

### Section 3: numbers presenting to hospital each day ###

# daily events
counts.per.day = group_by(transitions, simnum, day, From, To) %>%
  summarise(n=n()) %>%
  arrange(simnum, n)
# numbers presenting to hospital that are infected or not
hospital = filter(transitions, To %in% c('I2','HS','HE','HR')) %>%
  mutate(hospital = case_when(To == 'I2' ~ 'Infected', # serious infection presenting to hospital
                              To == 'H1' ~ 'Infected', # mild infection presenting to hospital
                              To == 'HE' ~ 'Infected', # exposed infection presenting to hospital
                              To == 'HS' ~ 'Not Infected', # susceptible and recovered
                              To == 'HR' ~ 'Not Infected')) %>%
  group_by(simnum, day, hospital) %>%
  summarise(n=n())

# compare daily ratio of infected to not infected
compare = tidyr::spread(hospital, hospital, n) %>%
  mutate(Infected = replace_na(Infected, 0),
         `Not Infected` = replace_na(`Not Infected`, 0),
         p = Infected / (Infected + `Not Infected`) ) %>%
  ungroup()
dsmooth = group_by(compare, day) %>%
  summarise(p = mean(p)) %>%
  mutate(simnum=51)
to.plot = bind_rows(compare, dsmooth) %>% select(simnum, day, p)
colours = c(grey(runif(n=50, min=0.2, max=0.8)), 'red') # line colours
lplot = ggplot(to.plot, aes(x=day, y=p, col=factor(simnum)))+
  geom_line(lwd=0.5)+
  scale_color_manual(NULL, values=colours)+
  xlab('Days')+
  ylab('P(infected)')+
  theme_bw()+
  coord_cartesian(xlim=c(0, 150))+
  theme(legend.position = 'none')
lplot 


## Section 4: age at transition ###

age.mean = group_by(transitions, simnum, From, To) %>%
  summarise(mean = mean(transitionAge), Q1=quantile(transitionAge, 0.25), Q3=quantile(transitionAge, 0.75))
# plot, first remove hospital transitions
for.plot = filter(age.mean,
                  !From %in% c('HR','HS','HE','H1','R','S'),
                  !To %in% c('HR','HS','HE','H1','S','dead')) %>% # remove some transitions to neaten plot
  mutate(
    FromNice = factor(From, levels=c('E','I1','I2','I3','R','D','HR','HS','H1','HE','dead'), # make nicer labels
                    labels = c('Exposed','Mild\ninfection','Severe\ninfection','Critical\ninfection','Recovered','Died - COVID','Hospital R','Hospital S','Hospital I1','Hospital E','Died - Other')),
    ToNice = factor(To, levels=c('E','I1','I2','I3','R','D','HR','HS','H1','HE','dead'), # make nicer labels
                         labels = c('Exposed','Mild infection','Severe infection','Critical infection','Recovered','Died - COVID','Hospital R','Hospital S','Hospital I1','Hospital E','Died - Other'))) 
age.box = ggplot(data=for.plot, aes(x=FromNice, y=mean, fill=ToNice))+
  geom_boxplot()+
  scale_fill_manual('Transition to', values=cbPalette)+
  ylab('Age')+
  xlab('Transition from')
age.box
#
jpeg('figures/mean.age.transition.jpg', width=5, height=4, units='in', res=300)
print(age.box)
dev.off()

