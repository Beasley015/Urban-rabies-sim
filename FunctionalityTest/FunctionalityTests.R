# Function tests for urban ABM

library(tidyverse)
library(ggnewscale)
library(patchwork)
library(gifski)
library(gganimate)
library(sf)

setwd("./FunctionalityTest")

# Movement/distance decay test -------------
raccoons <- read.csv("mvt_test.csv") %>%
  filter(rep == sample(1:max(rep), 1))

hr_props <- raccoons %>%
  filter(week >= 27 & week <= 40) %>%
  group_by(id, x, y, hab) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(prop = count/sum(count)) %>%
  filter(prop > 0.1) %>%
  ungroup()
  
hr_sizes <-hr_props %>%
  dplyr::select(id, x, y) %>%
  group_by(id, x, y) %>%
  distinct() %>%
  st_as_sf(coords = c("y", "x")) %>%
  ungroup() %>%
  group_by(id) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_buffer(dist = 0.5) %>%
  st_concave_hull(ratio = 0.3, allow_holes = F)
  
summary(st_area(hr_sizes))

# Look at land cover-------------
# Load in data
landscape <- read.csv(file='ExampleLand.csv', header = F)

# Convert to long format
df <- landscape %>%
  as.data.frame(.) %>%
  mutate(y = 1:nrow(landscape)) %>%
  pivot_longer(cols = -y, names_to = 'x', values_to = 'hab') %>%
  mutate(x = rep(1:ncol(landscape), nrow(landscape))) %>%
  mutate(hab = as.character(hab)) %>%
  mutate(hab = factor(hab, levels = c(1:7, 10)))

levels(df$hab) = c("Deciduous", "DevLo", "Pasture",
                                 "DevHi", "Wetlands", "Conifers",
                                 "Crops", "Buffer")

ggplot(data = df, aes(x=x, y=y, fill=hab))+
  geom_tile()+
  scale_fill_viridis_d(name = "Habitat", begin = 1, end = 0)+
  theme(axis.title = element_blank(), axis.text = element_blank(),
        legend.position = "None")
        #legend.key.size = unit(0.2, 'in'))

# ggsave(filename = "sample_landscape.jpeg", height = 2, width = 4,
#        units = "in")

# Length of latent period ------------------
dis <- read.csv("disease.csv")

ggplot(data = dis, aes(x = time_since_inf))+
  geom_bar(fill = 'lightgray', color = 'black')+
  labs(x = "Length of latent period (Weeks)", y = "Count")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

ggsave(filename = "latent_period.jpeg", height = 3, width = 4,
       units = "in")

# Recovery probability -----------
rec <- read.csv("disease.csv")

rec2 <- rec %>%
  group_by(week) %>%
  summarise(n_transition = n(), n_recover = min(recover)) %>%
  ungroup() %>%
  summarise(perc = sum(n_recover)/sum(n_transition))

rec2

# Mortality tests RESUME HERE --------------------
mort_counts <- read.csv(file = "mvt_mortality_test.csv")

# Mortality as % of population
perc_df <- mort_counts %>%
  mutate(year = c(rep(1, 52), rep(2,52))) %>%
  rowwise(step) %>%
  mutate(perc_deaths = sum(n_random_mort, n_dis_mort, orphan_mort,
                           juvie_cc_mort, adult_cc_mort)/
           pop_size)

summary(perc_df$perc_deaths)
summary(perc_df$perc_deaths[perc_df$year==1])
summary(perc_df$perc_deaths[perc_df$year==2])

ggplot(perc_df, aes(x = step, y = perc_deaths, 
                    color = factor(year)))+
  geom_point(size = 2)+
  scale_color_manual(values = c("black", "limegreen"),
                     name = "Year")+
  labs(x = "Week", y = "Mortality Rate")+
  theme_bw(base_size = 13)+
  theme(panel.grid = element_blank())

ggsave(filename="mortality_Test.jpeg", height = 3, width = 6,
       units = "in")

# Causes of mortality
more_percs <- mort_counts %>%
  mutate(year = c(rep(1, 52), rep(2,52))) %>%
  rowwise(step) %>%
  mutate(total_deaths = sum(n_random_mort, n_dis_mort, orphan_mort,
                            juvie_cc_mort, adult_cc_mort)) %>%
  mutate(across(n_random_mort:adult_cc_mort, 
                .fns = ~ .x/total_deaths,
                .names = "{.col}_perc"))

no_dis <- more_percs %>%
  filter(year == 1) %>%
  select(n_random_mort_perc:adult_cc_mort_perc) #%>%
  # filter(step > 20 & step < 40)

summary(no_dis$juvie_cc_mort_perc)
colMeans(no_dis, na.rm = T)

dis <- more_percs %>%
  filter(year == 2) %>%
  select(n_random_mort_perc:adult_cc_mort_perc)

colMeans(dis, na.rm = T)

