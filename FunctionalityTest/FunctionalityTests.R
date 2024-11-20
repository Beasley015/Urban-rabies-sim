# Function tests for urban ABM

library(tidyverse)
library(ggnewscale)
library(patchwork)

setwd("./FunctionalityTest")

# Test some distance-decay functions -----------
x <- seq(0,3000,100)
y <- exp(-0.0005 * (x^2)/100)

plot(x,y)

# Movement/distance decay test -------------
raccoons <- read.csv("mvt.csv")

guys.to.keep <- sample(raccoons$id, 5)

a.few.raccoons <- raccoons %>%
  filter(id %in% guys.to.keep, week < 43) %>%
  group_by(id, x, y, habs) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(prop = count/sum(count))

home.coords <- raccoons %>%
  filter(week == 1)

ggplot(data = a.few.raccoons, aes(x = x, y = y, fill = prop,
                                  color = factor(id)))+
  geom_tile(size = 1)+
  scale_fill_gradient(low = 'white', high = 'forestgreen',
                      name = "Proportion")+
  scale_color_viridis_d(name = "ID")+
  geom_point(data = home.coords[home.coords$id %in% guys.to.keep,],
             aes(x = x, y = y), size = 4, color = "black",
             inherit.aes = F)+
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
landscapes <- array(0, dim = c(60,60,3))

for(i in 1:3){
  tab <- as.matrix(read.csv(file= paste('example_land', 
                                               i, '.csv', sep = ""),
                              header = F))
  
  landscapes[,,i] <- tab
}

plotlist <- list()

for(i in 1:3){
# Convert to long format
  df <- landscapes[,,i] %>%
    as.data.frame(.) %>%
    rename_with(.fn = ~as.character(1:ncol(landscapes))) %>%
    mutate(y = 1:nrow(landscapes)) %>%
    pivot_longer(cols = -y, names_to = 'x', values_to = 'hab') %>%
    mutate(x = as.numeric(x)) %>%
    mutate(hab = factor(hab))

  plotlist[[i]] <- ggplot(data = df, aes(x=x, y=y, fill=hab))+
    geom_tile()+
    scale_fill_viridis_d(name = "Habitat", begin = 1, end = 0)+
    theme(axis.title = element_blank(), axis.text = element_blank(),
          legend.key.size = unit(0.2, 'in'))
}

plotlist[[1]] + plotlist[[2]] + plotlist[[3]] +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'a')

# ggsave(filename = "sample_landscapes.jpeg", height = 2, width = 6,
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
