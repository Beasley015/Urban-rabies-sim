library(tidyverse)
library(viridis)

setwd("./ParamSensitivity")

dat <- read.csv("carrying_capacity_mortality.csv") %>%
  select(rep, year, week, total_pop, ac_mort, jc_mort) %>%
  mutate(nweek = ((year-1)*52)+week)

ggplot(data=dat, aes(x = nweek, y = total_pop, color = factor(rep)))+
  geom_line()+
  facet_grid(rows = vars(ac_mort), cols=vars(jc_mort))+
  scale_color_viridis_d(end = 0.9)+
  # geom_vline(xintercept = ((0:9)*52)+20, linetype = "dashed")+
  # geom_vline(xintercept = ((0:9)*52)+43, linetype = "dotted")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("full_cc_birthrate95.jpeg", width=10, height=9, units = "in")

jmort_only <- dat %>%
  filter(ac_mort == 0.005) %>%
  group_by(jc_mort, nweek) %>%
  summarise(mean_pop = mean(total_pop))

ggplot(data=jmort_only, aes(x = nweek, y = mean_pop, 
                            color = factor(jc_mort)))+
  geom_line(size = 1.5)+
  scale_color_viridis_d(end = 0.9)+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("cc_adult005_birthrate95.jpeg", width=5, height=3, units = "in")
