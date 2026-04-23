# Script to Quantity Assessment
#---------------------------------
# Load phylogeny and dataset
#---------------------------------
library(ape); library(phylolm); library(dplyr); library(tidyr)
tree <- read.tree("data/Phylo_tree_EUNIS.tree")
data <- read.csv("data/Appendix S6.csv") %>%
  mutate(species = gsub(" ", "_", Accepted_binomial)) %>%
  separate(species, into = c("genus", "spp."), sep = "_") %>%
  mutate(species = paste(genus, spp., sep = "_")) %>%
  select(-Accepted_binomial) %>%
  filter(species %in% tree$tip.label)# Full dataset of the 1,495 species analysed, including trait data and quality assessment results 

# Drop tips not in dataset
tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))

#---------------------------------------------------------------------
# Function to Phylogenetic logistic model for categorical predictors
#---------------------------------------------------------------------
run_phyloglm <- function(data, variables, predictor_group_name, phylo_tree, response = "GermData_availability") {
  results <- list()
  models <- list()
  
  for (var in variables) {
    df <- data %>%
      filter(!is.na(.data[[response]]), !is.na(.data[[var]]), species %in% phylo_tree$tip.label)
    
    rownames(df) <- df$species
    
    tree_pruned <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, rownames(df)))
    if (length(tree_pruned$tip.label) < 3) next
    
    model <- phyloglm(as.formula(paste(response, "~", var)), data = df, phy = tree_pruned, method = "logistic_MPLE")
    
    results[[var]] <- summary(model)$coefficients %>%
      as.data.frame() %>%
      slice(-1) %>%
      mutate(Predictor = var, Alpha = model$alpha)
    
    models[[var]] <- model
  }
  
  summary_df <- bind_rows(results) %>%
    add_row(Predictor = predictor_group_name, .before = 1) %>%
    select(Predictor, everything())
  
  return(list(summary = summary_df, models = models))
}


#---------------------------------------------------------------------
# Function to Phylogenetic logistic model for continuous predictors
#---------------------------------------------------------------------
run_phyloglm_cont <- function(data, variables, predictor_group_name, phylo_tree, response = "GermData_availability", log_transform = TRUE) {
  results <- list()
  models <- list()
  
  for (var in variables) {
    df <- data %>%
      filter(!is.na(.data[[response]]), !is.na(.data[[var]]), species %in% phylo_tree$tip.label)
    
    rownames(df) <- df$species
    
    tree_pruned <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, rownames(df)))
    
    if (length(tree_pruned$tip.label) < 3) next
    
    formula <- if (log_transform) {
      as.formula(paste(response, "~ log(", var, ")"))
    } else {
      as.formula(paste(response, "~", var))
    }
    
    model <- phyloglm(formula, data = df, phy = tree_pruned, method = "logistic_MPLE")
    
    coef_df <- summary(model)$coefficients %>%
      as.data.frame() %>%
      slice(-1) %>%
      mutate(Predictor = var,
             Alpha = model$alpha)
    
    results[[var]] <- coef_df
    models[[var]] <- model
  }
  
  summary_df <- bind_rows(results) %>%
    add_row(Predictor = predictor_group_name, .before = 1) %>%
    select(Predictor, everything())
  
  return(list(summary = summary_df, models = models))
}

#----------------------------------------------
# Example 1: Run model for categorical variables: 
#----------------------------------------------
## *Lifeform 
dat_life <- data %>%
  mutate(perennial = ifelse(Lifeform == "perennial", 1, 0),
         woody = ifelse(Lifeform == "woody", 1, 0),
         annual = ifelse(Lifeform == "annual", 1, 0))%>%
  mutate(perennial=scale(perennial),
         woody = scale(woody),
         annual = scale(annual))

rownames(dat_life) <- dat_life$species

lifeform <- run_phyloglm(dat_life, c("perennial", "woody", "annual"), "Lifeform", tree)
lifeform_df <- lifeform$summary

## Pollination 
## Fruit type 
## Dormancy 
## Conservation interest
## Restoration relevance
## Reported uses

#----------------------------------------------
# Example 2: Run model for continuous variables
#----------------------------------------------
## *Seed mass
mass <- run_phyloglm_cont(data,variables = c("Seed_mass"),
                          predictor_group_name = "Seed mass",
                          phylo_tree = tree)
mass_df <- mass$summary

## Ecological breath



