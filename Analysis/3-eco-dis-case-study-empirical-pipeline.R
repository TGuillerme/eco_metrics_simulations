



############################################################################################################################
############################################ ECODISSIMILARITY CASE STUDY: HAWAII ###########################################
############################################################################################################################

setwd("C:/Users/mwj207/OneDrive - University of Birmingham/hawaii/eco-dissimilarity/github")

library(devtools)
require(treats)
require(dispRity)
require(BAT)
require(FD)
require(ks)
require(TPD)
library(tidyverse)

## Loading the functions
source("./Functions/dbFD.R") # I added this in order to import the dbFD function
source("./Functions/melodic.rao.R")
source("./Functions/make.simulation.pipeline.R")
source("./Functions/simulation.pipeline.R")

## Loading the optimisation functions
source("./Functions/optim.replicate.R")
source("./Functions/optim.helpers.R")
source("./Functions/analyse.replicates.R")



###### SLIGHTLY ALTERED FUNCTIONS ##########################################################################################

fun.proba.den.2 = function(presence, traits, verbose) {
  if(verbose) {
    message(".", appendLF = FALSE)
  }
  
  if(ncol(traits) <= 4) {
    traits_4dmax <- traits
  } else {
    traits_4dmax <- traits[, c(1:4)] # moved this from under "Measure probability density metrics
  }
  
  ## Get the sds matrix
  standard_deviations <- matrix(
    data = rep(sqrt(diag(ks::Hpi.diag(traits_4dmax))), nrow(presence)), # Hpo.diag was originally: traits
    nrow = nrow(presence),
    byrow = TRUE,
    dimnames = list(rownames(presence), paste0("sdA", 1:ncol(traits_4dmax)))) # ncol argument was origninally: traits
  
  ## Measure probability density metrics

  suppressMessages(TPD_species <- TPD::TPDsMean(species = rownames(presence),
                                                means = traits_4dmax,
                                                sds = standard_deviations,
                                                alpha = 0.99))
  TPD_cases   <- TPD::TPDc(TPDs = TPD_species, sampUnit = t(presence))
  output <- suppressMessages(do.call(cbind, TPD::REND(TPDc = TPD_cases)$communities))
  colnames(output) <- paste("TPD", gsub("F", "", colnames(output)), sep = "_")
  return(output)
}


simulation.pipeline.2 = function(sim.data, type, remove, verbose = FALSE, record.timer = FALSE, n.taxa) { # n.taxa argument it added to provide number of species in original assemblage
  
  if(record.timer) {
    time_start <- Sys.time()
  }
  
  if(verbose) {
    message("\nSimulating the trait space:", appendLF = FALSE)
  }
  
  if(is(sim.data, c("matrix")) || is(sim.data, c("data.frame"))) {
    ## Sim.data is an input traitspace
    is_empirical <- TRUE
    ## Check if remove is the correct list
    if(!is(remove, "list") || !is(unlist(remove), "logical")) {
      stop("If providing a empirical traitspace, you must provide a list of logicals for remove.")
    }
  } else {
    is_empirical <- FALSE
  }
  
  timer <- NULL
  
  ## Make empty reductions (all FALSE) for checking if the reduction worked
  if(!is_empirical) {
    target_reductions <- replicate(length(remove), FALSE, simplify = FALSE)
    time_start <- Sys.time()
    while(any(!unlist(lapply(target_reductions[which(!(remove == 0 & remove == 1))], check.reduction)))) {
      
      if(verbose) {
        message(".", appendLF = FALSE)
      }
      
      ## Simulating the brownian birth death data
      simulated_data <- treats::treats(
        bd.params = list(speciation = sim.data$speciation),
        stop.rule = list(max.living = sim.data$n.taxa),
        traits    = treats::make.traits(process = treats::BM.process, n = sim.data$n.traits))
      trait_space <- simulated_data$data[rownames(simulated_data$data) %in% simulated_data$tree$tip.label, ]
      
      if(verbose) {
        message(".", appendLF = FALSE)
      }
      
      ## Applying the stressors
      random_reductions <- sapply(remove, function(remove, trait_space) return(dispRity::reduce.space(trait_space, type = "random", remove = remove)), trait_space, simplify = FALSE)
      
      if(verbose) {
        message(".", appendLF = FALSE)
      }
      
      ## Run the reductions
      target_reductions <- switch(type,
                                  "facilitation" = sapply(remove, function(remove, trait_space) return(dispRity::reduce.space(trait_space, type = "density", remove = remove)), trait_space, simplify = FALSE),
                                  "equalizing"   =  sapply(remove, function(remove, trait_space) return(dispRity::reduce.space(trait_space, type = "size", remove = remove)), trait_space, simplify = FALSE),
                                  "competition"  =  sapply(remove, function(remove, trait_space) return(dispRity::reduce.space(trait_space, type = "evenness", remove = remove, parameters = list(power = 3))), trait_space, simplify = FALSE),
                                  "filtering"    =  sapply(remove, function(remove, trait_space) return(dispRity::reduce.space(trait_space, type = "position", remove = remove)), trait_space, simplify = FALSE)
      )
    }
    
    ## Wrapping up the data
    presences <- mapply(function(random, stressor) return(cbind(random, stressor)), random = random_reductions, stressor = target_reductions, SIMPLIFY = FALSE)
    
    if(record.timer) {
      time_end <- Sys.time()
      timer <- c("Data simulation" = time_end - time_start)
    }
    
    if(verbose) message("Done.\n", appendLF = FALSE)
    
  } else {
    ## Get trait_space
    trait_space <- sim.data
    if(is.null(rownames(trait_space))) {
      rownames(trait_space) <- 1:n.taxa # changed to 1: n.taxa to allow for changing the number of species in original assemblage 
    }
    ## Create the presences list
    make.presences <- function(one_remove) {
      ## Create the random variable (null)
      random <- rep(FALSE, length(one_remove))
      random[sample(1:n.taxa, size = sum(one_remove))] <- TRUE # changed to 1: n.taxa to allow for changing the number of species in original assemblage 
      ## Get the dimnames
      if(is.null(names(one_remove))) {
        rownames <- 1:length(one_remove)
      } else {
        rownames <- names(one_remove)
      }
      return(matrix(c(random, one_remove), ncol = 2, byrow = FALSE, dimnames = list(c(rownames), c("random", "stressor"))))
    }
    presences <- lapply(remove, make.presences)
  }
  
  ## Calculating the functional dendrograms
  if(verbose) {
    message("Calculating functional dendrogram:", appendLF = FALSE)
  }
  
  ## Making the trees  
  tree <- hclust(dist(trait_space), method = "average")
  
  ## Run the tree analyses
  tree_results <- lapply(presences, fun.tree, tree = tree, verbose)
  
  if(record.timer) {
    time_end <- Sys.time()
    timer <- c(timer, "Tree metrics" = time_end - time_start)
    time_start <- Sys.time()
  }
  
  ## Calculating the functional dissimilarity based metrics
  if(verbose) {
    message("Done.\n", appendLF = FALSE)
    message("Calculating dissimilarity:", appendLF = FALSE)
  }
  
  ## Run the FD analyses
  FD_results <- lapply(presences, fun.dissimilarity, traits = trait_space, verbose)
  
  if(record.timer) {
    time_end <- Sys.time()
    timer <- c(timer, "Dissimilarity metrics" = time_end - time_start)
    time_start <- Sys.time()
  }
  
  ## Calculating the convex hull
  if(verbose) {
    message("Done.\n", appendLF = FALSE)
    message("Calculating convex hull:", appendLF = FALSE)
  }
  
  ## Run the convex.hull analyses
  convex_hull_results <- lapply(presences, fun.convex.hull, traits = trait_space, verbose)
  
  if(record.timer) {
    time_end <- Sys.time()
    timer <- c(timer, "Convex hull metrics" = time_end - time_start)
    time_start <- Sys.time()
  }
  
  ## Calculating the kernel density hypervolumes
  if(verbose) {
    message("Done.\n", appendLF = FALSE)
    message("Calculating kernel density hypervolumes:", appendLF = FALSE)
  }
  
  ## Run the convex.hull analyses
  hypervolume_results <- lapply(presences, fun.kernel, traits = trait_space, verbose)
  # hypervolume_results <- lapply(presences, fun.tree, tree = tree, verbose) ; warning("DEBUG")
  
  if(record.timer) {
    time_end <- Sys.time()
    timer <- c(timer, "Hypervolume metrics" = time_end - time_start)
    time_start <- Sys.time()
  }
  
  ## Calculating the trait probabilities density
  if(verbose) {
    message("Done.\n", appendLF = FALSE)
    message("Calculating trait probability densities:", appendLF = FALSE)
  }
  
  ## Get the results
  TPD_results <- lapply(presences, fun.proba.den.2, traits = trait_space, verbose) # changed fun.propa.den.2 function to work for 5 axes 
  # TPD_results <- lapply(presences, fun.tree, tree = tree, verbose) ; warning("DEBUG")
  
  if(verbose) {
    message("Done.\n", appendLF = FALSE)
  }
  
  if(record.timer) {
    time_end <- Sys.time()
    timer <- c(timer, "Probability density metrics" = time_end - time_start)
  }
  
  ## Combine all the results
  all_results <- list(tree_results, FD_results, convex_hull_results, hypervolume_results, TPD_results)
  while(length(all_results) > 1) {
    all_results[[1]] <- mapply(cbind, all_results[[1]], all_results[[2]], SIMPLIFY = FALSE)
    all_results[[2]] <- NULL
  }
  all_results <- all_results[[1]]
  
  if(record.timer) {
    return(list("results" = all_results, "trait_space" = trait_space, "reductions" = presences, "timer" = timer))
  } else {
    return(list("results" = all_results, "trait_space" = trait_space, "reductions" = presences))
  }
}




make.simulation.pipeline.2 = function(type = NULL,
                                        sim.data = list(n.traits   = 5, # set to 5 dimensions
                                                        speciation = 1,
                                                        n.taxa     = n.taxa), # add argument for number of species 
                                        remove = c((63-37)/63), # set random reductions to the same reductions as for the empirical examples 
                                        verbose  = TRUE,
                                        record.timer = TRUE, 
                                      n.taxa = n.taxa) {
  return(function() simulation.pipeline.2(sim.data, type, remove, verbose = verbose, record.timer = record.timer, n.taxa = n.taxa))
}




###### LOAD TRAIT SPACE ####################################################################################################

## global trait space 
empirical_traitspace = read_csv('../trait-space-5d-prehistoric-community.csv') 

## historic species list 
hist.spec = read_csv('../historic-species-list-minus-a-flammeus.csv') %>%  
  select(species )

## extant species list 
pres.spec = read_csv('../present-species-list-minus-a-flammeus.csv') %>% 
  select(species) 


## prehistoric trait space 
prehist_traitspace = empirical_traitspace %>% 
  column_to_rownames('species') %>% 
  as.matrix() 

## historic trait space 
hist_traitspace = empirical_traitspace %>% 
  filter(empirical_traitspace$species %in% hist.spec$species)%>% 
  column_to_rownames('species') %>% 
  as.matrix() 




######## LIST OF SPECIES TO REMOVE FOR EACH SUBSAMPLE ######################################################################

# A list of species to remove for each sub-sampling the list should be list of named `logical` vectors 
# (where for example species that should be excluded - e.g. fossils are marked as `FALSE` and the ones to include as `TRUE`).

ph.subsamples = list(ifelse(rownames(prehist_traitspace) %in% hist.spec$species, TRUE, FALSE)) # P>H: (118-63)/118 = 47% reduction 
pc.subsamples = list(ifelse(rownames(prehist_traitspace) %in% pres.spec$species, TRUE, FALSE)) # P>C (118-38)/118 = 68% reduction
hc.subsamples = list(ifelse(rownames(hist_traitspace) %in% pres.spec, TRUE, FALSE)) # H>C(63-37)/63 = 40% 


prehist.subsamples = c(ph.subsamples, pc.subsamples) # P>H & P>C


## Adding the elements names
prehist.subsamples = lapply(prehist.subsamples, function(x) {names(x) <- rownames(prehist_traitspace) ; return(x)})
str(prehist.subsamples)

hist.subsamples = lapply(hc.subsamples, function(x) {names(x) <- rownames(hist_traitspace) ; return(x)})
str(hist.subsamples)




######## RUNNING THE EMPIRICAL PIPELINE ####################################################################################
## Creating the function for running the empirical pipeline
my.empirical.pipeline <- make.simulation.pipeline.2(sim.data = prehist_traitspace,
                                                    remove   = prehist.subsamples, 
                                                    n.taxa = 64
                                                    )

## Running the pipeline
my.empirical.pipeline()



## Running the full optimisation
results.ext <- optim.replicate(input.fun = my.empirical.pipeline,
                           diagnose = var,
                           summarise = optim.summarise,
                           minimum = 20, maximum = 170,
                           stop.variance = 0.0001,
                           verbose = TRUE,
                           bkp.path = "../eco-matrix-processed/",
                           bkp.name = "empirical_data-hist-extant.rda") 



