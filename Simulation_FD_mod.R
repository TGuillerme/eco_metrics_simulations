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
library("FD")
library("geometry")
library("hypervolume")
library("TPD")

## ------------------------------------------------------------------------
# Custom R functions ------------------------------------------------------
## ------------------------------------------------------------------------

## Function for simplifying the plots (1D)
#' @param data the dataset
#' @param reduction the vector of selected element for the reduction
#' @param null the vector of selected element for the null
#' @param col two colours
#' @param pch the pch (default = 19)
#' @param trait which trait to plot (default is 1 for plot.1D and c(1,2) for plot.2D)
#' @param xlab the xlab
#' @param legend the legend text
plot.1D <- function(data, reduction, null, main, col = c("blue", "orange"), pch = 19, trait = 1, xlab = "trait values", legend) {
  ## Get the histograms
  main_plot <- hist(data, plot = FALSE)
  plot_1 <- hist(data[null, trait], plot = FALSE, breaks = main_plot$breaks)
  plot_2 <- hist(data[reduction, trait], plot = FALSE, breaks = main_plot$breaks)
  plot_height <- max(plot_1$counts, plot_2$counts)
  
  ## Make the colours as transparent (for overlaps)
  make.transparent <- function(col, alpha = 0.5) {
    rgb_vals <- grDevices::col2rgb(col)
    rgb_args <- list(red   = rgb_vals["red", ],
                     green = rgb_vals["green", ],
                     blue  = rgb_vals["blue", ],
                     maxColorValue = 255,
                     alpha = alpha * 255)
    return(do.call(grDevices::rgb, rgb_args))
  }
  ##Plot the histograms
  plot(plot_1, col = make.transparent(col[1]), ylim = c(0, plot_height),
       main = ifelse(missing(main), "", main),
       xlab = "Trait values", ylab = "Counts")
  plot(plot_2, col = make.transparent(col[2]), add = TRUE)
  
  ## Adding the legend
  if(!missing(legend)) {
    legend(x = "topleft", col = col, pch = pch, legend = legend)
  }
  ## Add the points
  points(x = data[null, trait], y = rep(0, sum(null)), col = "black", bg = col[1], pch = 21)
  points(x = data[reduction, trait], y = rep(-plot_height/50, sum(reduction)), col = "black", bg = col[2], pch = 21)
}

## Function for simplifying the plots 2D
#' @param data the dataset
#' @param reduction the vector of selected element for the reduction
#' @param null the vector of selected element for the null
#' @param col two colours
#' @param pch the pch (default = 19)
#' @param trait which trait to plot (default is 1 for plot.1D and c(1,2) for plot.2D)
#' @param xlab the xlab
#' @param legend the legend text
plot.2D <- function(data, reduction, null, main, col = c("blue", "orange"), pch= 19, trait = c(1,2), xlab = "trait values", legend) {
  ## Plotting the histogram
  plot(NULL, xlim = range(trait_space[, trait[1]]), ylim = range(trait_space[, trait[2]]), main = main)
  points(trait_space[null, trait], pch = pch, col = col[1])
  points(trait_space[reduction, trait], pch = pch, col = col[2])
  ## Adding the legend
  if(!missing(legend)) {
    legend(x = "topleft", col = col, pch = pch, legend = legend)
  }
}

## ------------------------------------------------------------------------
# Simulating the pool trait space -----------------------------------------
## ------------------------------------------------------------------------

# We will build simulate data using a Brownian Motion (the breath of trait values just increases with time, unconstrained), 
# a pure birth tree (no extinction) and will stop the simulation as soon as we reach 200 taxa. 
# We stick to 2 dimensions (2 traits) to ease visualization

## Creating the trait object (2D BM)
my_trait <- make.traits(process = BM.process, n = 2)

## Setting up the birth death parameters (no extinction)
bd_params <- list(speciation = 1)

## Setting up the stopping rules
stop_rule <- list(max.living = 200)

## Simulating the data
simulated_data <- dads(bd.params = bd_params,
                       stop.rule = stop_rule,
                       traits    = my_trait)

## We can then visualise the simulated data (just to show the process):
par(mfrow = c(2, 1))
plot(simulated_data, main = "The simulated BM through time")
plot(simulated_data, main = "The 2D BM", trait = c(1,2), xlab = "traiy1", ylab = "trait2")

## The orange data points are node values and will be ignored
## Here we'll select only the tip values from the data
trait_space <- simulated_data$data[rownames(simulated_data$data) %in% simulated_data$tree$tip.label, ]

## ------------------------------------------------------------------------
# Simulating the different mechanisms -------------------------------------
## ------------------------------------------------------------------------

# The stressors here are the processes we apply to change the normal pattern of the data
# (the BM motion generated normal distribution of the data) 
# for simulating the different mechanisms. 
# To simulate these patterns, we can use the function dispRity::reduce.space that has already several algorithms implemented. 

# Null mechanism ----------------------------------------------------------
#This is the really basic one: just a null model were we randomly sample 100 species out of the 200. 
#This can be done using the stats::sample function but for making the whole script neater, we can also use the "random" algorithm in dispRity::reduce.space that just randomly removes n% of species:

## Removing 50% of the elements randomly
random_reduction <- reduce.space(trait_space, type = "random", remove = 0.5)

# And we can visualise these results compared to the full distribution of 200 species:
par(mfrow = c(2, 1))
plot.1D(trait_space, reduction = random_reduction, null = rep(TRUE, 200), main = "1D null mechanism", legend = c("Full distribution (200sp)", "Random reduction (100sp)"))
plot.2D(trait_space, reduction = random_reduction, null = !random_reduction, main = "2D null mechanism")

# Environmental filtering -------------------------------------------------
# For the environmental filtering algorithm, we can use the We can use the "position" algorithm from the dispRity::reduce.space 
# function that shifts the position of the group to one corner of the space so that X% of the data in on side and (1-X)% is on the other side 
# (with side being weirder and weirder as dimensionality increases!):
   
## Removing 50% to shift them
filtering_reduction <- reduce.space(trait_space, type = "position", remove = 0.5)

## Visualising the results
par(mfrow = c(2, 1))
plot.1D(trait_space, reduction = filtering_reduction, null = random_reduction, main = "1D filtering mechanism")
plot.2D(trait_space, reduction = filtering_reduction, null = random_reduction, main = "2D filtering mechanism")

# Competitive exclusion ---------------------------------------------------
# For simulating the competitive exclusion, 
# we use the dispRity::reduce.space "evenness" algorithm that basically flattenshe curve (i.e displacing things from the highest density regions to the sides).

## Removing 50% to flatten the curve
competitive_exclusion <- reduce.space(trait_space, type = "evenness", remove = 0.5, parameters = list(power = 3))

## Visualising the results
par(mfrow = c(2, 1))
plot.1D(trait_space, reduction = !competitive_exclusion, null = random_reduction, main = "1D excluding mechanism")
plot.2D(trait_space, reduction = !competitive_exclusion, null = random_reduction, main = "2D excluding mechanism")

# Equalising fitness processes --------------------------------------------
#For simulating equalisation of fitness we can use the "size" algorithm from dispRity::reduce.space. 
# This removes elements around past a distance from the centre of the population 
# (with the distance being determined automatically so that X% are inside the radius and 1-X% are outside):
  
## Removing 50% to change the trait space size
equalizing_fitness <- reduce.space(trait_space, type = "size",remove = 0.5)

## Visualising the results
par(mfrow = c(2, 1))
plot.1D(trait_space, reduction = equalizing_fitness, null = random_reduction, main = "1D equalising mechanism")
plot.2D(trait_space, reduction = equalizing_fitness, null = random_reduction, main = "2D equalising mechanism")

# Facilitation ------------------------------------------------------------
# We use the "density" algorithm from dispRity::reduce.space that will reduce the nearest neighbours distance between elements 
# (basically making clumps of two elements close to each other).

## Removing 50% to change the trait space size
facilitation <- reduce.space(trait_space, type = "density", remove = 0.5)

## Visualising the results
par(mfrow = c(2, 1))
plot.1D(trait_space, reduction = facilitation, null = random_reduction, main = "1D facilitation mechanism")
plot.2D(trait_space, reduction = facilitation, null = random_reduction, main = "2D facilitation mechanism")

## ------------------------------------------------------------------------
# Testing the performance of FD methods -----------------------------------
## ------------------------------------------------------------------------

# Now that we have our mechanisms, we can run simulation to test the performance of different functional diversity methods in
# detecting change in the trait space (see Mammola et al. 2020 EcoEvoRxiv, doi: 10.32942/osf.io/j64nt).

# We run two set of analysis: 
# 1) comparing the properties of the trait space generated by the mechanism (richness, divergence and regularity) against the full trait space.
# 2) comparing the properties of the trait space generated by the mechanism (richness, divergence and regularity) against 
#    a randomly sampled trait space (i.e. null modelling)

# Reference trait space ---------------------------------------------------

# Storing everything we simulated before
data <- data.frame(trait_space,
                   full = rep(TRUE,200),
                   random_reduction,
                   filtering_reduction,
                   competitive_exclusion,
                   equalizing_fitness,
                   facilitation)

#converting true/false in 0/1
data[,3:8] <- ifelse(data[,3:8] == TRUE,1,0)

# Dissimilarity based metrics
dbFD <- FD::dbFD(trait_space, print.pco = FALSE, calc.FGR = FALSE) 

FD_Rao        <- dbFD$RaoQ
FD_divergence <- dbFD$FDiv
FD_evenness   <- dbFD$FEve

# Functional dendrogram
tree <- hclust(dist(trait_space))
tree_richness   <- BAT::alpha(comm = rep(1,200), tree = tree)
tree_dispersion <- BAT::dispersion(comm = rep(1,200), tree = tree) 
tree_regularity <- BAT::evenness(comm = rep(1,200), tree = tree)

par(mfrow = c(1, 1)) ; plot(tree)

# Convex hull
hull <- geometry:::convhulln(trait_space, options = "FA")
hull_richness <- BAT::hull.alpha(comm = hull)
# dispersion and evenness not possible (see Mammola et al. 2020 EcoEvoRxiv, doi: 10.32942/osf.io/j64nt)

coords <- trait_space[ c( chull(trait_space) , chull(trait_space) [1]), ]
par(mfrow = c(1, 1)) ; plot(trait_space, xlab = "Trait1", ylab = "Trait2") ; lines(coords, col = "blue") #visualizing it

# Kernel density hypervolume
hv <- hypervolume::hypervolume_gaussian(trait_space)
hv_richness   <- BAT::kernel.alpha(comm = hv)
hv_dispersion <- BAT::kernel.dispersion(comm = hv)
hv_regularity <- BAT::kernel.evenness(comm = hv)  

par(mfrow = c(1, 1)) ; plot(hv, col="blue") #visualizing it

# Trait Probability Density
bw <- sqrt(diag(ks::Hpi.diag(data[, 1:2])))
spSD <- matrix(rep(bw, nrow(data)), nrow = nrow(data), byrow = T, 
               dimnames = list(rownames(data), c("sdA1", "sdA2")))
TPDSpecies <- TPDsMean(species = rownames(data),
                       means = data[, 1:2],
                       sds = spSD,
                       alpha = 0.99)
TPDcases <- TPDc(TPDs = TPDSpecies, sampUnit = t(data[,3:8]))

#Richness, evenness and divergence of each case:
TPDResults <- REND(TPDc = TPDcases)$communities

TPD_richness   <- TPDResults$FRichness[1]
TPD_dispersion <- TPDResults$FDivergence[1]
TPD_regularity <- TPDResults$FEvenness[1]

par(mfrow = c(1, 1)) ; plotTPD(TPDcases,whichPlot=1) # visualizing it

# Starting the simulation -------------------------------------------------

n.perm <- 5 #Number of permutations 

#empty data.frame
result <-  data.frame(simulation   = NULL,
                      method_class = NULL,
                      method       = NULL,
                      mechanism    = NULL,
                      Richness     = NULL,
                      Dispersion   = NULL,
                      Evenness     = NULL)

for (p in 1 : n.perm) {
  
  message(paste("Starting simulation", p, "out of", n.perm, sep = ' '))
  
  # Simulating the different processes
  random_reduction      <- reduce.space(trait_space, type = "random", remove = 0.5)
  facilitation          <- reduce.space(trait_space, type = "density", remove = 0.5)
  equalizing_fitness    <- reduce.space(trait_space, type = "size", remove = 0.5)
  competitive_exclusion <- reduce.space(trait_space, type = "evenness", remove = 0.5, parameters = list(power = 3))
  filtering_reduction   <- reduce.space(trait_space, type = "position", remove = 0.5)
  
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
  
  tree_p <- hclust(dist(trait_space_p))
  tree_richness   <- BAT::alpha(comm = t(mechanism_p), tree = tree_p)
  tree_dispersion <- BAT::dispersion(comm = t(mechanism_p), tree = tree_p) 
  tree_regularity <- BAT::evenness(comm = t(mechanism_p), tree = tree_p)

  tree_result <- cbind(tree_richness,tree_dispersion,tree_regularity)
  
  # Estimating Dissimilarity based metrics
  message("------ Estimating Dissimilarity based metrics ------")
  
  dbFD <- FD::dbFD(x = trait_space_p, a = t(mechanism_p), print.pco = FALSE, calc.FGR = FALSE) 
  
  FD_Rao        <- dbFD$RaoQ
  FD_divergence <- dbFD$FDiv
  FD_evenness   <- dbFD$FEve
  
  Rao_result  <- cbind(rep(NA,5),FD_Rao,rep(NA,5))
  Fdiv_result <- cbind(rep(NA,5),FD_divergence,rep(NA,5))
  FD_evenness <- cbind(rep(NA,5),rep(NA,5),FD_evenness)
  
  colnames(Rao_result) <- colnames(Fdiv_result) <- colnames(FD_evenness) <- colnames(tree_result)
  
  # Convex hull
  message("------ Estimating convex hull ------")
  
  hull_richness <- BAT::hull.alpha(comm = t(mechanism_p), trait = trait_space_p)
  # dispersion and evenness not possible (see Mammola et al. 2020 EcoEvoRxiv, doi: 10.32942/osf.io/j64nt)
  
  hull_result <- cbind(hull_richness,rep(NA,5),rep(NA,5))
  colnames(hull_result) <- colnames(tree_result)
  rownames(hull_result) <- rownames(tree_result)
  
  # Kernel density hypervolume
  message("------ Estimating kernel density hypervolumes ------")
  
  hv <- BAT::kernel.alpha(comm = t(mechanism_p), trait = trait_space_p, return.hv = TRUE)
  hv_richness   <- hv[[1]]
  hv_dispersion <- BAT::kernel.dispersion(comm = hv[[2]])
  hv_regularity <- BAT::kernel.evenness(comm = hv[[2]])  
 
  hv_result <- t(rbind(hv_richness,hv_dispersion,hv_regularity))
  colnames(hv_result) <- colnames(tree_result)
  
  # Trait Probability Density
  message("------ Estimating trait probability densities ------")
  
  # bw <- sqrt(diag(ks::Hpi.diag(trait_space_p)))
  # spSD <- matrix(rep(bw, nrow(mechanism_p)), nrow = nrow(mechanism_p), byrow = T, 
  #                dimnames = list(rownames(mechanism_p), c("sdA1", "sdA2")))
  # TPDSpecies <- TPDsMean(species = rownames(mechanism_p),
  #                        means = trait_space_p,
  #                        sds = spSD,
  #                        alpha = 0.99)
  # TPDcases <- TPDc(TPDs = TPDSpecies, sampUnit = t(mechanism_p))
  
  # #Richness, evenness and divergence of each case:
  # TPDResults <- REND(TPDc = TPDcases)$communities
   
  # #Storing the result
  # TPD_result <- t(rbind(TPDResults$FRichness,TPDResults$FDivergence,TPDResults$FEvenness))
  # colnames(TPD_result) <- colnames(tree_result)
  
  ## Storing all tresult of the permutation
  result_p <- rbind(Rao_result,Fdiv_result,FD_evenness,tree_result,hull_result,hv_result)#,TPD_result)
  
   result2 <-  data.frame(simulation   = rep(paste("Sim_", p, sep = ''), nrow(result_p)),
                      method_class = c(rep("dissimilarity",15),rep("tree",5),rep("binary hypervolume",5),
                                       rep("probabilistic hypervolume",5)),
                      method       = c(rep("Rao",5),rep("Fdiv",5),rep("FEve",5),
                                       rep("tree",5),rep("hull",5),rep("kernel",5)),
                      mechanism    = rownames(result_p),
                      result_p)
  
   result <- rbind(result,result2)
} 



