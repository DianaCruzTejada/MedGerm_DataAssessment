# Thermal modelling (Sentinella et al 2020)
library(glmmTMB); library(tidyverse); library(MCMCglmm); library(writexl)
read.csv("data/original.csv") %>% # from Cruz-Tejada et al., 2024b
  select(id_test,accepted_binomial,Family, clade5,seedlot,Tmean,Germinated,Germinable) -> data
length(unique(data$accepted_binomial))# 459 spp

#Prepare data ####
data <- data %>% # dataset name 
  mutate(propgerm = Germinated/Germinable )%>%   #new column with the germination proportion
  group_by(accepted_binomial) %>%
  mutate(maxgerm = max(propgerm)) %>%
  data.frame()
#Find maximum germination percentage in group
data <- data %>%
  group_by(accepted_binomial) %>%
  mutate(maxgerm = max(propgerm))

###### Number of tested temperatures (NTestTemp > 2) #### 
# Flagged any groups with less than 3 different temperatures as 0
pf1 <- data %>%
  group_by(accepted_binomial) %>%
  add_tally(n_distinct(Tmean), name = "Tmean_count")%>%
  mutate(pf1=ifelse(Tmean_count > 2,1,0))%>%
  select(accepted_binomial,pf1)%>%
  distinct()

# filter out the species that do not have more than two data of temperature
MCMC <- data %>%
  group_by(accepted_binomial) %>%
  add_tally(n_distinct(Tmean), name = "Tmean_count")%>%
  filter(any(Tmean_count > 2))

length(unique(MCMC$accepted_binomial))#322

# Modeling ####
safe_glm <- safely(glmmTMB)

TMB_models <- MCMC %>% group_by(accepted_binomial) %>%
  do(germ.model = safe_glm(cbind(Germinated, Germinable - Germinated)~ I(Tmean) + I((Tmean) ^ 2)
                           + (1|seedlot)+(1|id_test),
                           ziformula = ~0,
                           dispformula=~1,family="betabinomial",data=.))

# load thermal modelling results
#load("results/TMB_models.RData")
#### Extract coefficients a, b, and c from models (parabola) ####
TMB_MaxMinOpt <- TMB_models %>%
  rowwise() %>%
  mutate(
    a = ifelse(is.null(germ.model[[1]]), NA, fixef (germ.model[[1]])$cond[3]),
    b = ifelse(is.null(germ.model[[1]]), NA, fixef(germ.model[[1]])$cond[2]),
    c = ifelse(is.null(germ.model[[1]]), NA, fixef(germ.model[[1]])$cond[1])
  )

##### Add in Topt, Tmax, Tmin estimates from coefficients ####
TMB_MaxMinOpt <- data.frame(
  TMB_MaxMinOpt %>%
    mutate(
      Topt = (-b / (2 * a)),
      Logit.Max = a * Topt * Topt + b * Topt + c,
      #Maximum germination from model at Topt (logit)
      Model.Max = (exp(Logit.Max) / (1 + exp(Logit.Max))),
      #convert to actual maximum germination
      Topt.upp = (-b - (sqrt((
        b * b - 4 * a * (c - log(0.95 * Model.Max / (1 - 0.95 * Model.Max)))
      )))) / (2 * a),
      Topt.low = (-b + (sqrt((
        b * b - 4 * a * (c - log(0.95 * Model.Max / (1 - 0.95 * Model.Max)))
      )))) / (2 * a),
      Tmax = (-b - (sqrt((
        b * b - 4 * a * (c - log(0.05 / (1 - 0.05)))
      )))) / (2 * a),
      Tmin = (-b + (sqrt((
        b * b - 4 * a * (c - log(0.05 / (1 - 0.05)))
      )))) / (2 * a),
      Tbreadth = Tmax - Tmin,
      Toptbreadth = Topt.upp - Topt.low
    ))


# Delta method - standard error estimate (functions) ####
se.matrix.fun <-
  function(a,
           b,
           c,
           Model.Max,
           prop,
           upp.low,
           germ.model,
           is.opt) {
    if (is.null(germ.model[[1]])) {NA}
    else {
      x <-
        ifelse(is.opt == T,
               log(prop * Model.Max / (1 - prop * Model.Max)), #If calulating Topt, use model maximum
               log(prop / (1 - prop))) #If calculating Min/Max, use total proportion
      delta <- b * b - 4 * c * a + 4 * x * a
      #Original error calculations
      pa <- ifelse(upp.low == "upp",
                   1 / (2 * a ^ 2) * (b + (delta ^ (1 / 2)) + 2 * a * (c - x) * (delta ^ (-1 / 2))), #Topt.upp/Tmax
                   1 / (2 * a ^ 2) * (b - (delta ^ (1 / 2)) - 2 * a * (c - x) * (delta ^ (-1 / 2)))) #Topt.low/Tmin
      pb <- ifelse(upp.low == "upp",
                   -1 / (2 * a) * (1 + b * (delta ^ (-1 / 2))), #Topt.upp/Tmax
                   -1 / (2 * a) * (1 - b * (delta ^ (-1 / 2))))#Topt.low/Tmin
      pc <- ifelse(upp.low == "upp", (delta ^ (-1 / 2)), #Topt.upp/Tmax
                   -(delta ^ (-1 / 2))) #Topt.low/Tmin
      vector.v <- as.vector(c(pc, pb, pa))
      vector.h <- t(vector.v) #transposed vector
      vcov.matrix <- as.matrix(vcov(germ.model[[1]])$cond)
      variance <- vector.h %*% vcov.matrix %*% vector.v
      variance <- ifelse(variance < 0, NA, sqrt(as.numeric(variance)))
    }
  }

se.opt.fun <- function(a, b, germ.model) {
  if (is.null(germ.model[[1]])) {
    NA
  }
  else {
    pa <- b / (2 * a * a)
    pb <- -1 / (2 * a)
    pc <- 0
    vector.v <- as.vector(c(pc, pb, pa))
    vector.h <- t(vector.v) #transposed vector
    vcov.matrix <- as.matrix(vcov(germ.model[[1]])$cond)
    vcov.matrix[1, ] <- 0
    vcov.matrix[, 1] <- 0
    variance <- vector.h %*% vcov.matrix %*% vector.v
    variance <- ifelse(variance < 0, NA, sqrt(as.numeric(variance)))
  }
}
se.breadth.fun <- function(a, b, c, Model.Max, prop, germ.model, is.opt) {
  if (is.null(germ.model[[1]])) {
    NA
  }
  else {
    x <-
      ifelse(is.opt == T,
             log(prop * Model.Max / (1 - prop * Model.Max)), #If calulating Topt breadth, use model maximum
             log(prop / (1 - prop))) #If calculating breadth, use total proportion
    delta <- b * b - 4 * c * a + 4 * x * a
    ## Corrected formula
    pa <- (-1 / (a * a)) * (sqrt(delta) + (2 * a * (c - x) / (sqrt(delta))))
    pb <- b / (a * sqrt(delta))
    pc <- -2 / sqrt(delta)
    vector.v <- as.vector(c(pc, pb, pa))
    vector.h <- t(vector.v) #transposed vector
    vcov.matrix <- as.matrix(vcov(germ.model[[1]])$cond)
    variance <- vector.h %*% vcov.matrix %*% vector.v
    variance <- ifelse(variance < 0, NA, sqrt(as.numeric(variance)))
  }
}
# Standard error #####
TMB_MaxMinOpt <- TMB_MaxMinOpt %>%
  rowwise() %>%
  mutate(
    SEmax = se.matrix.fun(a, b, c, Model.Max, 0.05, "upp", germ.model, F),
    SEmin = se.matrix.fun(a, b, c, Model.Max, 0.05, "low", germ.model, F),
    SEopt.upp = se.matrix.fun(a, b, c, Model.Max, 0.95, "upp", germ.model, T),
    SEopt.low = se.matrix.fun(a, b, c, Model.Max, 0.95, "low", germ.model, T),
    SEopt = se.opt.fun(a, b, germ.model),
    SEbreadth = se.breadth.fun(a, b, c, Model.Max, 0.05, germ.model, F),
    SEoptbreadth = se.breadth.fun(a, b, c, Model.Max, 0.95, germ.model, T)
  )

# Remove models ####
TMB_MaxMinOpt <- dplyr::select(TMB_MaxMinOpt,-germ.model)
res_data <- left_join(TMB_MaxMinOpt, MCMC, "accepted_binomial")
res_data <- as.data.frame(res_data)

## FLAGGING ####
#### Table 2. Thermal modelling criteria used to assess the germination data quality.
res_data <- merge(pf1,res_data,by="accepted_binomial",all.x = T)

## Model concavity: a < 0 & !is.na(a) ####
res_data <- res_data %>%
  mutate(Af1 = if_else(a < 0 & !is.na(a), 1, 0))

## Model flatness: TestTemp > Topt & propgerm < 0.75*Pmax | TestTemp< Topt & propgerm < 0.75*Pmax
res_data <- res_data %>% group_by(accepted_binomial) %>%
  mutate(Af7 = ifelse((any((Tmean > Topt & propgerm < 0.75 * maxgerm)) &
                         any((Tmean < Topt & propgerm < 0.75 * maxgerm))), 1, 0))

# remove duplicates
res_data <- res_data %>%
  select(-id_test,-seedlot,-Tmean,-Germinated,-Germinable,-propgerm,-maxgerm)%>%
  distinct()

##### Extract results ####
# write.csv(res_data, "results/res_ThermalModelling_version_139.csv")
write.csv(res_data, "results/res_ThermalModeling.csv")





