library("FactoMineR"); library("factoextra"); library(tidyverse);library(grid);library(gridExtra);library(ggimage)

# ==================================
# Load data
# ==================================
data <- read.csv("results/Supporting Elements/Appendix_S6.csv") # This dataset contains thermal modelling results 

# ==================================
## Thermal modeling criteria: species passed quality assesment 
# ==================================
d1 <- data %>%
  filter(GermData_quality == 1)

# ==================================  
## Mean of Topt and Toptbreath of Mediterranean species
# ==================================
Tmean_opt <- round(mean(d1$Topt)) #17°C
Tmean_opt.breath <- round(mean(d1$Toptbreadth)) #8°C

# ==================================
## Risk allocation
# ==================================
#Thermally resilient --> Topt > 17°C and Toptbreath > 8°C
#Warm-narrow risk --> Topt > 17°C and Toptbreath < 8°C
#Cold-wide risk --> Topt < 17°C and Toptbreath > 8°C
#Thermally costrained --> Topt < 17°C and Toptbreath < 8°C

d2 <- d1 %>%
  mutate(risk = ifelse(Topt >17 & Toptbreadth>8,"Thermally resilient",
                       ifelse(Topt > 17& Toptbreadth <=8,"Warm-narrow",
                              ifelse(Topt <= 17 & Toptbreadth > 8,"Cold-wide",
                                     ifelse(Topt <= 17 & Toptbreadth <= 8,"High", "NADA")))))


