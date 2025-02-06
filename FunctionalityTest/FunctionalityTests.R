# Function tests for urban ABM

library(tidyverse)
library(ggnewscale)
library(patchwork)
library(gifski)
library(gganimate)
library(sf)

setwd("./FunctionalityTest")

# Test some distance-decay functions -----------
x <- seq(0,3000,100)
y <- exp(-0.001 * (x^2)/100)

plot(x,y)

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

guys.to.keep <- sample(raccoons$id[which(raccoons$week==1)], 5)

a.few.raccoons <- raccoons %>%
  filter(id %in% guys.to.keep, week < 43) %>%
  group_by(rep, id, x, y, hab) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(prop = count/sum(count))

home.coords <- raccoons %>%
  filter(week == 1)

few.raccoons.sp <- dplyr::select(a.few.raccoons, id, x, y) %>%
  filter(id %in% c(348, 62))
coordinates(few.raccoons.sp) <- c("x", "y")

raccoon.mcp <- mcp(few.raccoons.sp, percent = 95)

ggplot(data = a.few.raccoons, aes(x = x, y = y, fill = prop,
                                  color = factor(id)))+
  geom_tile(size = 1)+
  # geom_polygon(data = raccoon.mcp, aes(long,lat, group = id))
  scale_fill_gradient(low = 'white', high = 'forestgreen',
                      name = "Proportion")+
  scale_color_viridis_d(name = "ID")+
  # geom_point(data = home.coords[home.coords$id %in% guys.to.keep,],
  #            aes(x = x, y = y), size = 4, color = "black",
  #            inherit.aes = F)+
  lims(x = c(0,21), y = c(0,21))+
  theme_dark(base_size = 14)

# ggsave(filename = "samplemvt.jpeg", width = 6, height=5,
#        units = "in")

# Look at effects of habitat -----------
# HR size based on habitat in attractor cell
hab.vals <- data.frame(id = 1:5, vals = c("Deciduous", "DevLo",
                                          "Pasture", "DevHi",
                                          "Wetlands"))
habitats <- raccoons %>%
  filter(week < 43) %>%
  group_by(id, x, y, habs) %>%
  summarise(count = n()) %>%
  left_join(hab.vals, by = c("habs" = "id"))

ggplot(data = habitats, aes(x = vals, y = count))+
  geom_boxplot(fill='lightgray')+
  labs(x = "Habitat", y = "Annual Home Range Size (# cells)")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave("sizebyhab.jpeg", width = 4, height = 3, units = "in")

summary(aov(count~factor(habs), data = habitats))

# Proportion of raccoons w/attractor per habitat type
hr_habs <- home.coords %>%
  group_by(habs) %>%
  summarise(count = n()) %>%
  mutate(Used = count/sum(count)) %>%
  mutate(Available = c(0.2585, 0.2337, 0.1915, 0.1266, 
                       0.0899)) %>%
  mutate(names = c("Deciduous", "DevLo", "Pasture", "DevHi", 
                   "Wetlands")) %>%
  select(-c(count, habs)) %>%
  pivot_longer(cols = -names, names_to = "var", values_to = "prop")

ggplot(data = hr_habs, aes(x = names, y = prop, fill = var))+
  geom_col(position = "dodge", color = 'black')+
  scale_fill_manual(values = c('limegreen', 'lightgray'), name = "")+
  labs(y = "Proportion", x = "Habitat")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave("HomeRangeHabs.jpeg", height = 4, width = 6)

# Time spent in each habitat type
hab_total <- raccoons %>%
  filter(week < 43) %>%
  group_by(hab_current) %>%
  summarise(count = n()) %>%
  mutate(Used = count/sum(count)) %>%
  mutate(Available = c(0.2585, 0.2337, 0.1915, 0.1266, 
                       0.0899)) %>%
  mutate(names = c("Deciduous", "DevLo", "Pasture", "DevHi", 
                   "Wetlands")) %>%
  select(-c(count, hab_current)) %>%
  pivot_longer(cols = -names, names_to = "var", values_to = "prop")

ggplot(data = hab_total, aes(x = names, y = prop, fill = var))+
  geom_col(position="dodge", color = 'black')+
  
  scale_fill_manual(values = c('limegreen', 'lightgray'), name = "")+
  labs(y = "Proportion", x = "Habitat")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave("PositionHabs.jpeg", height = 4, width = 6)

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

ggplot(data = dis, aes(x = time))+
  geom_bar(fill = 'lightgray', color = 'black')+
  labs(x = "Length of latent period (Weeks)", y = "Count")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

ggsave(filename = "latent_period.jpeg", height = 3, width = 4,
       units = "in")

# Recovery probability -----------
rec <- read.csv("recovery.csv")

rec2 <- apply(rec, 2, sum)/sum(rec)

# Disease spread gif and figs ------------
# Set up data
dis <- read.csv("mvt_disease_test.csv") %>%
  mutate(nweek = ((year-1)*52)+ week) %>%
  group_by(nweek, x, y) %>%
  summarise(inc = n(), inf = sum(inf)) %>%
  mutate(inf = case_when(inf == 0 ~ '0',
                         TRUE ~ '1'))

# Create gif
trans.fig <- ggplot(data = dis, aes(x = x, y = y, fill = inc))+
  geom_tile(aes(width = 1, height = 1), linewidth = 1)+
  geom_point(aes(x = x, y = y, alpha = inf), size = 2,
             color = 'white')+
  scale_fill_viridis_c(name = "Exposed")+
  scale_alpha_discrete(range = c(0, 1), name = "Infectious")+
  lims(x = c(0, 21), y = c(0, 21))+
  theme(axis.title = element_blank(), 
        title = element_text(size = 16))

animation_test <- 
  trans.fig + 
  transition_states(nweek)+
  labs(subtitle = "Week: {closest_state}") 

trans.animate <- animate(animation_test, nframes = 45, fps = 2)
trans.animate

# anim_save(filename="trans_animation.gif", 
#           animation = trans.animate)

# Create stationary figures
dis_subset <- filter(dis, nweek <= min(nweek)+3)
sub.figs <- list()

for(i in 1:length(unique(dis_subset$nweek))){
  sub.figs[[i]] <- ggplot(data = dis_subset[dis_subset$nweek==unique(dis_subset$nweek)[i],], aes(x = x, y = y, fill = inc))+
    geom_tile(aes(width = 1, height = 1), linewidth = 1)+
    geom_point(aes(x = x, y = y, alpha = inf), size = 2,
               color = 'white')+
    scale_fill_viridis_c(name = "Infected", limits = c(0,max(dis_subset$inc)))+
    scale_alpha_discrete(range = c(0, 1), name = "Infectious",
                         guide = "none")+
    lims(x = c(0, 21), y = c(0, 21))+
    labs(title = paste("Week", unique(dis_subset$nweek)[i], 
                       sep = " "))+
    theme_bw()
    theme(axis.title = element_blank(), 
          title = element_text(size = 18))
}

(sub.figs[[1]] + sub.figs[[2]] + sub.figs[[3]])+
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = 'a')

# ggsave(filename = "fixed_trans_fig.jpeg", width = 12, height = 4,
#        units = "in")

# Mortality tests --------------------
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
