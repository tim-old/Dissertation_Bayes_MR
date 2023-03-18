library(tidyverse)
library(janitor)

alcohol_data <- read.csv("alcohol_consumption.csv")

# reformat data 
alcohol_data_tidy <- alcohol_data %>% 
  pivot_wider(names_from = c(Variable, Measure),
              values_from = Value,
              id_cols = c(Year, Country))


glimpse(alcohol_data_tidy)
skimr::skim(alcohol_data_tidy) 
names(alcohol_data_tidy)
