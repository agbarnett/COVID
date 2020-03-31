# 0_aus_pop.R
# get the Australian population data from the ABS
# data from https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/3101.0Dec%202018?OpenDocument
# March 2020
library(readxl)
library(dplyr)
infile = read_excel('data/3101059.xls', sheet='Data1', skip=10, col_names = FALSE) %>%
  mutate(`...1` = as.character(`...1`)) %>%
  filter(`...1` == '2018-06-01') %>%
  select(-`...1`) %>%
  tidyr::gather(value='count') %>%
  mutate(index = 1:n())
# extract men
males = filter(infile, index <= 101) %>%
  mutate(age = 1:n() - 1,
         sex = 'm')
# extract women
females = filter(infile, index > 101, index<=202) %>%
  mutate(age = 1:n() - 1,
         sex = 'f')
# combine
ages = bind_rows(males, females) %>%
  select(-index, -key) %>%
  filter(age>=18) %>% # just adults
  mutate(p = count / sum(count)) # probability for each yearly age and sex
# save
save(ages, file='data/AusPop.RData')
         