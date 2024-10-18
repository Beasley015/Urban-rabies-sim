# Function tests for urban ABM

library(tidyverse)
library(ggnewscale)

setwd("./FunctionalityTest")

# Movement/distance decay test -------------
raccoons <- read.csv("mvt.csv")

guys.to.keep <- sample(raccoons$id, 2)

a.few.raccoons <- raccoons %>%
  filter(id %in% guys.to.keep) %>%
  group_by(id, x, y) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(prop = count/sum(count))

home.coords <- raccoons %>%
  filter(week == 1)

ggplot()+
  geom_tile(data = a.few.raccoons[a.few.raccoons$id == 
                                    guys.to.keep[1],],
            aes(x = x, y = y, fill = prop))+
  scale_fill_gradient(as.character(guys.to.keep[1]), low = "white",
                      high = "red")+
  new_scale_fill()+
  geom_tile(data = a.few.raccoons[a.few.raccoons$id ==
                                    guys.to.keep[2],],
            aes(x = x, y = y, fill = prop))+
  scale_fill_gradient(as.character(guys.to.keep[2]), low = "white",
                      high = "purple")+
  geom_point(data = home.coords[home.coords$id %in% guys.to.keep,],
             aes(x = x, y = y), size = 4, color = "black",
             inherit.aes = F)+
  lims(x = c(0,20), y = c(0,20))+
  theme_dark(base_size = 14)

# ggsave(filename = "samplemvt.jpeg", width = 6, height=5, 
#        units = "in")
