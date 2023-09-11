# Earthworm abundance and diversity
# EES 4th year field course 2023
# Jess, Caitlin, Simon, Cory, Heather, Imara, Alex
# 11-12 Sept 2023 

# libraries ----
library(skimr)
library(tidyverse)
install.packages("scatterplot3d")
library(scatterplot3d)
library(vegan)
library(devtools)
devtools::install_github("cmartin/ggConvexHull")

# loading data ----
worms <- read.csv("worms/earthworm_data.csv")

# viewing and cleaning data ----
view(worms)
skim(worms)

worms_sub <- filter(worms, habitat %in% c("unmowed_grassland","mowed_grassland"))  # removing woodland data
view(worms_sub)

# raw data visualization ----
(rough <- ggplot(worms_sub, aes(x = habitat, y = earthworm_count)) +
  geom_boxplot() +
  theme_classic())  # any difference between worm counts in 2 habitats

(rough2 <- ggplot(worms_sub, aes(x = soil_ph, y = earthworm_count, color = habitat)) +
    geom_point() +
    theme_classic())  # pH impact on earthworm count

# Data wrangling ----

(worms_sub <- worms_sub %>% 
  mutate(epigeic_prop = epigeic / adults, 
         anecic_prop = anecic / adults,
         endogeic_prop = endogeic / adults))  # adding columns for ecotype proportions

(worms_sub <- worms_sub %>% 
    mutate(adult_prop = adults / earthworm_count, 
           juv_prop = juveniles / earthworm_count))  # adding columns for age proportions


# modelling ----

## modelling worm abundance as impacted by habitat
a_mod <- glm(earthworm_count ~ habitat, data = worms_sub, family = poisson)
summary(a_mod)
anova <- aov(a_mod)
plot(a_mod)


### CHECK FOR OVERDISPERSION

## modelling worm communities by proportions
### checking proportion distributions
hist(worms_sub$epigeic_prop, breaks = 30)
hist(worms_sub$anecic_prop, breaks = 30)
hist(worms_sub$endogeic_prop, breaks = 30)  # looks normal enough that we'll run with it


# plotting community structure by UG vs MG ----

## setting up data
comm <- worms_sub %>% 
  mutate(nonepi_prop = anecic_prop + endogeic_prop, 
         nonepi = anecic + endogeic) %>% 
  select(-c(date, soil_t1, soil_t2, soil_t3, soil_t4, soil_m1, soil_m2, soil_m3,
            soil_m4, anecic_prop, endogeic_prop, adult_prop, juv_prop))

## by counts
(community <- ggplot(comm, aes(x = epigeic, y = nonepi, col = habitat)) +
  geom_point() +
  ggConvexHull::geom_convexhull(alpha = 0.3, aes(fill = habitat)) +
  theme_classic())

## by proportions
(community <- ggplot(comm, aes(x = epigeic_prop, y = nonepi_prop, col = habitat)) +
    geom_point() +
    geom_polygon() +
    theme_classic())

# bar chart for community structure

## changing data to long format
comm_long <- comm %>% 
  pivot_longer(cols = c('epigeic', 'endogeic', 'anecic'),
               names_to = "ecotype", 
               values_to = "count")
View(comm_long)

## creating bar chart
(community_bar <- ggplot(comm_long, aes(x = habitat, y = count, fill = ecotype)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic())

# saving figures ----
ggsave("worms/img/community_struc_bar.png", plot = community_bar, width = 5, height = 5)
