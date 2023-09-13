# Earthworm abundance and diversity
# EES 4th year field course 2023
# Jess Thomson, Caitlin Carmichael, Simon Schulze, Cory MacCormack-Montequin, 
# Heather Young, Imara Thorpe, Alex Colety
# 11-12 Sept 2023 

# libraries ----
library(skimr)
library(tidyverse)
library(devtools)
devtools::install_github("cmartin/ggConvexHull")
library(ggConvexHull)

# loading data ----
worms <- read.csv("data/earthworm_data.csv")

# viewing and cleaning data ----
view(worms)
skim(worms)

worms_sub <- filter(worms, habitat %in% c("unmowed_grassland",
                                          "mowed_grassland"))  # removing woodland data
view(worms_sub)

# raw data visualization ----
(rough <- ggplot(worms_sub, aes(x = habitat, y = earthworm_count)) +
  geom_boxplot() +
  theme_classic())  # any difference between worm counts in 2 habitats?

(rough2 <- ggplot(worms_sub, aes(x = soil_ph, y = earthworm_count, color = habitat)) +
    geom_point() +
    theme_classic())  # pH impact on earthworm count?

# Data wrangling ----

(worms_sub <- worms_sub %>% 
  mutate(epigeic_prop = epigeic / adults, 
         anecic_prop = anecic / adults,
         endogeic_prop = endogeic / adults,
         adult_prop = adults / earthworm_count, 
         juv_prop = juveniles / earthworm_count))  # adding columns for ecotype 
                                                   # and age proportions

(comm <- worms_sub %>% 
  mutate(nonepi_prop = anecic_prop + endogeic_prop, 
         nonepi = anecic + endogeic) %>% 
  select(-c(date, soil_t1, soil_t2, soil_t3, soil_t4, soil_m1, soil_m2, soil_m3,
            soil_m4, anecic_prop, endogeic_prop, adult_prop, juv_prop)))  
  # creating data subset with anecic and endogeic worms grouped

(comm_long <- comm %>% 
  pivot_longer(cols = c('epigeic', 'endogeic', 'anecic', 'juveniles'),
               names_to = "ecotype", 
               values_to = "count"))  # changing data to long format


## checking distributions of ecotypes
### checking proportion distributions
hist(worms_sub$epigeic_prop, breaks = 30)
hist(worms_sub$anecic_prop, breaks = 30)
hist(worms_sub$endogeic_prop, breaks = 30)  # looks normal enough 


# modelling ----

## modelling worm abundance as impacted by habitat

abund_mod <- glm(earthworm_count ~ habitat, data = worms_sub, family = poisson)
summary(a_mod)
plot(a_mod)  # significant difference between worm abundance in 2 habitats, due 
             # to juvenile abundance differences

juv_hab_mod <- glm(juveniles ~ habitat, data = worms_sub, family = poisson)
summary(a_mod)
plot(a_mod)  # checking just abundance for juveniles

epi_hab_model <- glm(epigeic ~ habitat, data = worms_sub)
summary(epi_hab_model)

non_epi_hab <- glm(nonepi ~ habitat, data = comm, family = poisson)
summary(non_epi_hab)
plot(non_epi_hab)


## modelling environmental factors against habitats
MoistL_model <- glm(soil_m_avg ~ habitat, data = worms_sub)
summary(MoistL_model)  # no sig diff in moisture bt habitats

Temphdel <- lm(soil_t_avg ~ habitat, data = worms_sub)
summary(TempL_model) # no sig diff bt soil temps bt habitats

Soil_ph_hmodel <- lm(soil_ph ~ habitat, data = worms_sub)
summary(Soil_ph_hmodel) # soil pH not sig diff bt habitats


## modelling environmental factors against worms
Tempecmodel <- glm(soil_t_avg ~ earthworm_count, data = worms_sub)
summary(Tempecmodel) # No sig diff bt worm abundances with

TempWADmodel <- glm(soil_t_avg ~ adults, data = worms_sub)
summary(TempWADmodel) # Relationship bt soil temp and adult abundance, 
                      # prob due to epigeic relationship

TempWJmodel <- glm(soil_t_avg ~ juveniles, data = worms_sub)
summary(TempWJmodel) # no relationship bt juveniles and temp

Soil_ph_ecmodel <- glm(soil_ph ~ earthworm_count, data = worms_sub)
summary(Soil_ph_ecmodel) # soil pH unrelated to general earthworm abund



# plotting community structure by UG vs MG ----

## plotting ecotype community structures by count
(community <- ggplot(comm, aes(x = epigeic, y = nonepi, colour = habitat))     +
    geom_point()                                                               + 
    ggConvexHull::geom_convexhull(alpha = 0.3, aes(fill = habitat))            +
    theme_classic())                                                           +
  xlab("Average number of adult epigeic worms per plot")                       +
  ylab("Average number of adult non-epigeic worms per plot")                   +
  scale_fill_manual(name = "Habitat", values=c("#8BB26C", "#E16B7B"), 
                    labels = c("Mowed grassland", "Unmowed grassland"))        +
  scale_colour_manual(name = "Habitat", values=c("#8BB26C", "#E16B7B"), 
                      labels = c("Mowed grassland", "Unmowed grassland"))


# Bar chart for community structure by ecotype and abundances ----

### creating smaller dataframe w standard errors
comm_long_errors <- comm_long %>% 
  group_by(habitat, ecotype) %>% 
  summarise(mean = mean(count), se = (sd(count)/sqrt(length(count)))) %>% 
  ungroup()

### bar chart with new data frame
(community_bar_err <- ggplot(comm_long_errors, aes(x = habitat, y = mean, 
                                                   fill = ecotype))            +
    geom_bar(stat = "identity", position = "dodge")                            +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se, y = mean), 
                  width=0.25, position = position_dodge(0.9))                  +
    theme_classic())                                                           +
    xlab("Habitat")                                                              + 
    ylab("Average number of worms per plot")                                     +
    scale_fill_manual(name = "Ecotype", values = c("#8BB26C", "#225235", 
                                                   "#E16B7B", "#A9BFC4"), 
                    labels = c("Anecic", "Endogeic", "Epigeic", "Juveniles"))  +
    scale_x_discrete(labels = c("Mowed grassland", "Unmowed grassland"))         +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

## boxplot to check
(community_box <- ggplot(comm_long, aes(x = habitat, y = count, 
                                        fill = ecotype))                       +
    geom_boxplot()                                                             +
    theme_classic())

# Plotting habitat worm abundance against soil pH ----

(habitat_ph <- ggplot(worms, aes(x = habitat, y = soil_ph, color = habitat))   +
    geom_boxplot() +
    scale_colour_manual(name = "Habitat", values=c("#225235", "#E16B7B", 
                                                   "#8BB26C"), 
                       labels = c("Mixed Woodland", "Mowed grassland",
                                  "Unmowed grassland"))                        +
    labs(x = "Habitat", y = "Soil pH") +
    scale_x_discrete(labels = c("Mixed woodland", "Mowed grassland", 
                               "Unmowed grassland"))                           +
    theme_classic())


# saving figures ----
ggsave("img/community_struc_bar.png", plot = community_bar_err, width = 7,
       height = 5)
ggsave("img/community_struc_scatter.png", plot = community, width = 7, height = 5)
ggsave("img/habitat_ph.png", plot = habitat_ph, width = 7, height = 5)
