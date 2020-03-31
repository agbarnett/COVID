# 2_combine_simulations.R
# combine the multiple simulations
# March 2020
library(dplyr)

all_results = NULL
for (k in 1:50){
  infile = paste('Z:/COVID/data/simresults.', k, '.RData', sep='') # from lyra
  load(infile)
  transitions = mutate(transitions, simnum=k) %>%
    select(-birthDate) %>% # not needed
    mutate(#transitionTime = as.Date(transitionTime), # to avoid warning
           transitionTime = as.numeric(transitionTime))
  all_results = bind_rows(all_results, transitions) # concatenate results
}
# remove handful of transitions from D
all_results = filter(all_results, From != 'D')

# save big file
save(all_results, meta, file='data/all_results_vary.RData')
