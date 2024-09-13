library(tidyverse)

setwd("./ParamSensitivity")

dat <- read.csv("carrying_capacity_mortality.csv") %>%
  select(rep, year, week, total_pop, ac_mort, jc_mort) %>%
  mutate(nweek = ((year-1)*52)+week)

ggplot(data=dat, aes(x = nweek, y = total_pop, color = factor(rep)))+
  geom_line()+
  facet_grid(rows = vars(ac_mort), cols=vars(jc_mort))
