## ------------------------------------------------------------------------
## 'PAPER TITLE'
## ------------------------------------------------------------------------

# Authors::
# Thomas Guillerme (guillert@tcd.ie)
# Stefano Mammola (stefanomammola@gmail.com)
# Carlos Carmona (perezcarmonacarlos@gmail.com)
# Pedro Cardoso (pedro.cardoso@helsinki.fi)

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.0.3) and R studio (v. 1.4.1103)

## ------------------------------------------------------------------------
# Loading R packages ------------------------------------------------------
## ------------------------------------------------------------------------

library("BAT")
library("devtools")
# devtools::install_github("TGuillerme/dads")
# devtools::install_github("TGuillerme/dispRity")
library("dads")     # needs to be version >= 0.1.4
## Loading required package: dispRity
library("dispRity") # needs to be version >= 1.5.7
library("doSNOW")
library("FD")
library("geometry")
library("hypervolume")
library("parallel")
library("tcltk")
library("TPD")
source("melodic.rao.R")


## Simulation variables
## Data simulations
n_traits <- 2
speciation_rate <- 1
n_taxa <- 200

## Simulating the data
simulated_data <- dads::dads(
  bd.params = list(speciation = speciation_rate),
  stop.rule = list(max.living = n_taxa),
  traits    = dads::make.traits(process = dads::BM.process, n = n_traits))
trait_space <- simulated_data$data[rownames(simulated_data$data) %in% simulated_data$tree$tip.label, ]



#TODO:

# 1- Check variance > 0.05 or n_runs < 20
# 2- Simulate data
# 3- Run random + one other reduction at 0.2, 0.4, 0.6, 0.8
# 4- Measure metrics (store as lists)
# 5- Write outputs as csv (or keep in R)






# Starting the simulation -------------------------------------------------

n.perm <- 5 #Number of permutations 
#TG: should be minimum 20 and then increase by 5% until overall variance does not change by less than 5%


#Setting the cores
n.par <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

cl <- snow::makeSOCKcluster(n.par-2)
snow::registerDoSNOW(cl)

  ## This for a progress bar
  pb <- tcltk::tkProgressBar(max = n.perm)
  progress <- function(n) tcltk::setTkProgressBar(pb, n)
  opts <- list(progress=progress)

#setting the remove parameter for the simualation

remove_par = 0.5 #I waa thinking to run each at a time
# from 0.1 to 0.9

#Starting the parallelized simulation

## Transform into a while loop per stressor


result <- foreach(p=1:n.perm,
                  .combine = 'rbind',
                  .options.snow=opts,
                  .packages=c("dispRity","FD","TPD","BAT")) %dopar% {  #%do% to debug


  

  # Simulating the different processes
  random_reduction      <- dispRity::reduce.space(trait_space, type = "random", remove = remove_par)
  facilitation          <- dispRity::reduce.space(trait_space, type = "density", remove = remove_par)
  equalizing_fitness    <- dispRity::reduce.space(trait_space, type = "size", remove = remove_par)
  competitive_exclusion <- dispRity::reduce.space(trait_space, type = "evenness", remove = remove_par, parameters = list(power = 3))
  filtering_reduction   <- dispRity::reduce.space(trait_space, type = "position", remove = remove_par)
  
  # Storing the data
  data_p <- data.frame(trait_space,
                       random_reduction,
                       filtering_reduction,
                       competitive_exclusion,
                       equalizing_fitness,
                       facilitation)
  
  #converting true/false in 0/1
  data_p[,3:7] <- ifelse(data_p[,3:7] == TRUE,1,0)
  
  #selecting only columns with sp present in at least one mechanism
  data_p <- data_p[(rowSums(data_p[,3:7])>0)==TRUE,] 

  trait_space_p <- data_p[,1:2]
  mechanism_p   <- data_p[,3:7]
  






  # Functional dendrogram
  message("------ Estimating Functional dendrogram ------")

  ## Analyses  
  tree_p <- hclust(dist(trait_space_p), method = 'average')
  tree_richness   <- BAT::alpha(comm = t(mechanism_p), tree = tree_p)
  tree_dispersion <- BAT::dispersion(comm = t(mechanism_p), tree = tree_p) 
  tree_regularity <- BAT::evenness(comm = t(mechanism_p), tree = tree_p)
  ## Output
  tree_result <- cbind(tree_richness, tree_dispersion, tree_regularity)
  
  # Estimating Dissimilarity based metrics
  message("------ Estimating Dissimilarity based metrics ------")
  
  ## Analyses  
  dbFD <- FD::dbFD(x = trait_space_p, a = t(mechanism_p), print.pco = FALSE, calc.FGR = FALSE) 
  FD_Rao        <- dbFD$RaoQ
  FD_divergence <- dbFD$FDiv
  FD_evenness   <- dbFD$FEve
  ## Output
  Rao_result  <- cbind(rep(NA,5),FD_Rao,rep(NA,5))
  Fdiv_result <- cbind(rep(NA,5),FD_divergence,rep(NA,5))
  FD_evenness <- cbind(rep(NA,5),rep(NA,5),FD_evenness)

  #~~~~~~~~~~~
  # Why are there NA's here?
  #~~~~~~~~~~~

  colnames(Rao_result) <- colnames(Fdiv_result) <- colnames(FD_evenness) <- colnames(tree_result)
  
  # Convex hull
  message("------ Estimating convex hull ------")  

  ## Analyses
  hull_richness <- BAT::hull.alpha(comm = t(mechanism_p), trait = trait_space_p)
  ## Output
  hull_result <- cbind(hull_richness,rep(NA,5),rep(NA,5))
  colnames(hull_result) <- colnames(tree_result)
  rownames(hull_result) <- rownames(tree_result)

  #~~~~~~~~~~~
  # Why are there NA's here?
  #~~~~~~~~~~~  

  # Kernel density hypervolume
  message("------ Estimating kernel density hypervolumes ------")
  
  ## Analyses
  hv <- BAT::kernel.alpha(comm = t(mechanism_p), trait = trait_space_p, return.hv = TRUE)
  hv_richness   <- hv[[1]]
  hv_dispersion <- BAT::kernel.dispersion(comm = hv[[2]])
  hv_regularity <- BAT::kernel.evenness(comm = hv[[2]])  
  ## Output
  hv_result <- t(rbind(hv_richness,hv_dispersion,hv_regularity))
  colnames(hv_result) <- colnames(tree_result)
  
  # Trait Probability Density
  message("------ Estimating trait probability densities ------")

  ## Analyses  
  bw <- sqrt(diag(ks::Hpi.diag(trait_space_p)))
  spSD <- matrix(rep(bw, nrow(mechanism_p)), nrow = nrow(mechanism_p), byrow = T, 
                 dimnames = list(rownames(mechanism_p), c("sdA1", "sdA2")))
  TPDSpecies <- TPD::TPDsMean(species = rownames(mechanism_p),
                         means = trait_space_p,
                         sds = spSD,
                         alpha = 0.99)
  TPDcases <- TPD::TPDc(TPDs = TPDSpecies, sampUnit = t(mechanism_p))
  TPDResults <- TPD::REND(TPDc = TPDcases)$communities
   
  ## Output
  TPD_result <- t(rbind(TPDResults$FRichness,TPDResults$FDivergence,TPDResults$FEvenness))
  colnames(TPD_result) <- colnames(tree_result)
  

  ## Combining all the results
  result_p <- rbind(Rao_result,Fdiv_result,FD_evenness,tree_result,hull_result,hv_result,TPD_result)
  
  result2 <-  data.frame(simulation   = rep(paste("Sim_", p, sep = ''), nrow(result_p)),
                      method_class = c(rep("dissimilarity",15),rep("tree",5),rep("binary hypervolume",5),
                                       rep("probabilistic hypervolume",10)),
                      method       = c(rep("Rao",5),rep("Fdiv",5),rep("FEve",5),
                                       rep("tree",5),rep("hull",5),rep("kernel",5),rep("TPD",5)),
                      mechanism    = rownames(result_p),
                      result_p)
  
   return(result2)
} 

