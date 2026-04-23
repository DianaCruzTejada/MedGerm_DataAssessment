# Script to Quality Assessment

#---------------------------------
# Load phylogeny and dataset
#---------------------------------
library(phylolm); library(tidyverse); library(ape); library(tibble); 

tree <- read.tree("data/Phylo_tree_EUNIS.tree") 

data_qual <- read.csv("data/Appendix S6.csv") %>%
  mutate(species = gsub(" ", "_", Accepted_binomial)) %>%
  separate(species, into = c("genus", "spp."), sep = "_") %>%
  mutate(species = paste(genus, spp., sep = "_")) %>%
  select(-Accepted_binomial) %>%
  filter(species %in% tree$tip.label) 

tree <- drop.tip(tree, setdiff(tree$tip.label, data_qual$species))

#---------------------------------------------------------------------
# Function to Phylogenetic logistic model for categorical predictors
#---------------------------------------------------------------------
run_phyloglm <- function(data, variables, predictor_group_name, phylo_tree, response = "GermData_quality") {
  results <- list()
  models <- list()
  
  for (var in variables) {
    df <- data %>%
      filter(!is.na(.data[[response]]), !is.na(.data[[var]]), species %in% phylo_tree$tip.label)
    rownames(df) <- df$species
    tree_pruned <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, df$species))
    
    if (length(tree_pruned$tip.label) < 3) {
      warning(paste("Not enough species for binary predictor:", var))
      next
    }
    
    formula <- as.formula(paste(response, "~", var))
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

#---------------------------------------------------------------------
# Function to Phylogenetic logistic model for continuous predictors
#---------------------------------------------------------------------
run_phyloglm_cont <- function(data, variables, predictor_group_name, phylo_tree, response = "GermData_quality", log_transform = TRUE) {
  results <- list()
  models <- list()
  
  for (var in variables) {
    df <- data %>%
      filter(!is.na(.data[[response]]), !is.na(.data[[var]]), species %in% phylo_tree$tip.label) %>%
      column_to_rownames("species")
    
    tree_pruned <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, rownames(df)))
    if (length(tree_pruned$tip.label) < 3) next
    
    fml <- as.formula(paste(response, "~", if (log_transform) paste0("log(", var, ")") else var))
    model <- phyloglm(fml, data = df, phy = tree_pruned, method = "logistic_MPLE")
    
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

#----------------------------------------------
# Example 1: Run model for categorical variables
#----------------------------------------------
## *Lifeform 
dat_life <- data_qual %>%
  mutate(perennial = ifelse(Lifeform == "perennial", 1, 0),
         woody = ifelse(Lifeform == "woody", 1, 0),
         annual = ifelse(Lifeform == "annual", 1, 0)) %>%
  mutate(across(c(perennial, woody, annual), scale))
rownames(dat_life) <- dat_life$species

# Run model
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
mass <- run_phyloglm_cont(data_qual,variables = c("Seed_mass"),
                          predictor_group_name = "Seed mass",
                          phylo_tree = tree)
mass_df <- mass$summary

## Ecological breath 