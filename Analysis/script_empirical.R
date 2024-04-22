## Libraries
require(treats)
require(dispRity)
require(BAT)
require(FD)
require(ks)
require(TPD)

## Loading the functions
source("../Functions/melodic.rao.R")
source("../Functions/make.simulation.pipeline.R")
source("../Functions/simulation.pipeline.R")
source("../Functions/optim.replicate.R")
source("../Functions/optim.helpers.R")
source("../Functions/analyse.replicates.R")
source("../Functions/dbFD.R")

## global trait space 
empirical_traitspace <- read.csv('../Data/Raw/trait-space-5d-prehistoric-community.csv', row.names = 1) 

## Species lists
prehistoric_species <- rownames(read.csv("../Data/Raw/prehistoric-species-list.csv"))
historic_species <- rownames(read.csv('../Data/Raw/historic-species-list-minus-a-flammeus.csv', row.names = 1))
extant_species <- rownames(read.csv('../Data/Raw/present-species-list-minus-a-flammeus.csv', row.names = 1))

## Traitspace groups
extant_group <- rownames(empirical_traitspace) %in% extant_species ## 37 species present
historic_group <- rownames(empirical_traitspace) %in% historic_species ## 63 species present
prehistoric_group <- !(rownames(empirical_traitspace) %in% c(extant_species, historic_species)) ## 55 species present
names(extant_group) <- names(historic_group) <- names(prehistoric_group) <- rownames(empirical_traitspace)

groups <- list("extant" = extant_group, "historic" = historic_group, "prehistoric" = prehistoric_group)

## Set simulation type
simulate <- make.simulation.pipeline(sim.data = empirical_traitspace, remove = groups)

## Running the full optimisation
results <- optim.replicate(input.fun = simulate,
                           diagnose = var,
                           summarise = optim.summarise,
                           minimum = 20, maximum = 170,
                           stop.variance = 0.0001,
                           verbose = TRUE,
                           bkp.path = "../Data/Processed/",
                           bkp.name = "empirical_data.rda")