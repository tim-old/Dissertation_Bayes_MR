library(tidyverse)

alcohol_data <- read.csv("etoh_consumption.csv")


# Factorise years and country
alcohol_data_tidy <- alcohol_data %>% 
  mutate(Year = factor(Year)) %>% 
  mutate(Country = factor(Country)) %>% 
  arrange(Country, Year) %>% 
  select(Country, Year, Value) #%>%
  
# move countries as a variable 
#  pivot_wider(names_from = Country, values_from = Value)



