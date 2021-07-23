#' @title Simulation pipeline
#'
#' @description One simulation pipeline
#'
#' @param sim.data A named list of parameters for the data simulation (default is \code{list{n.traits = 2, speciation = 1, n.taxa = 200}})
#' @param type Which removal type (either \code{"facilitation"}, \code{"equalizing"}, \code{"competition"} or \code{"filtering"})}) to be paired with \code{"random"}.
#' @param remove How much data to remove (can be a vector of multiple percentages)
#' @param verbose whether to be verbose (\code{TRUE}, default) or not (\code{FALSE})
#' @param record.timer whether to record the timer for each step of the simulations
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

# sim.data = list(n.traits = 2, speciation = 1, n.taxa = 200)
# type = "facilitation"
# remove = c(0.2,.4, 0.6 ,0.8)
# library(dads)
# library(dispRity)
# library(BAT)
# source("melodic.rao.R")


# test <- simulation.pipeline(sim.data = sim.data, type = "facilitation", remove = remove, verbose = TRUE, record.timer = TRUE)

## SOLVED: Errors: Error in hclust(as.dist(thisTree)) : must have n >= 2 objects to cluster
# set.seed(4)
# test <- simulation.pipeline(sim.data, type = "facilitation", remove = remove, verbose = TRUE)
simulation.pipeline <- function(sim.data, type, remove, verbose = FALSE, record.timer = FALSE) {

    if(record.timer) {
        time_start <- Sys.time()
    }

    if(verbose) {
        message("\nSimulating the trait space:", appendLF = FALSE)
    }

    ## Make empty reductions (all FALSE) for checking if the reduction worked
    target_reductions <- replicate(length(remove), FALSE, simplify = FALSE)
    while(any(!unlist(lapply(target_reductions[which(!(remove == 0 & remove == 1))], check.reduction)))) {

        if(verbose) {
            message(".", appendLF = FALSE)
        }

        ## Simulating the brownian birth death data
        simulated_data <- dads::dads(
                bd.params = list(speciation = sim.data$speciation),
                stop.rule = list(max.living = sim.data$n.taxa),
                traits    = dads::make.traits(process = dads::BM.process, n = sim.data$n.traits))
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
        time_start <- Sys.time()
    }

    ## Calculating the functional dendrograms
    if(verbose) {
        message("Done.\n", appendLF = FALSE)
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
    TPD_results <- lapply(presences, fun.proba.den, traits = trait_space, verbose)
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

## Function for checking the reduction (needs at least 3 TRUE AND 3 FALSE)
check.reduction <- function(reduction) {
    return(length(which(reduction)) > 3 && length(which(!reduction)) > 3)
}
    
## Calculating the alpha, dispersion and eveness
fun.tree <- function(comm, tree, verbose) {
        if(verbose) {
            message(".", appendLF = FALSE)
        }
        table <- list(
            "tree_alpha"      = BAT::alpha(t(comm), tree),
            "tree_dispersion" = BAT::dispersion(t(comm), tree),
            "tree_evenness"   = BAT::evenness(t(comm), tree))
        table <- do.call(cbind, table)
        colnames(table) <- c("tree_alpha", "tree_dispersion", "tree_evenness")
        return(table)
}

## Calculating all the dissimilarity metrics
fun.dissimilarity <- function(presence, traits, verbose) {
    if(verbose) {
        message(".", appendLF = FALSE)
    }

    ## Only keep columns data have data in common
    keep_rows <- rowSums(presence) > 0
    ## Get the distance based metrics
    distance_FD <- FD::dbFD(x = traits[keep_rows, ], a = t(presence[keep_rows, ]), print.pco = FALSE, calc.FGR = FALSE, messages = FALSE)

    ## Get the melodic rao
    melo_rao <- melodic.rao(samp = t(presence[keep_rows, ]), dis = as.matrix(FD::gowdis(traits[keep_rows, ])))

    return(cbind("FD_Rao"        = distance_FD$RaoQ,
                 "FD_divergence" = distance_FD$FDiv,
                 "FD_evenness"   = distance_FD$FEve,
                 "melodic_Rao"   = melo_rao$rao,
                 "melodic_mpd"   = melo_rao$mpd))
}

## Calculating the hypervolume metrics
fun.kernel <- function(presence, traits, verbose) {
    if(verbose) {
        message(".", appendLF = FALSE)
    }

    ## Get the hypervolume (remove messages)
    suppressMessages(silent <- capture.output(hypervolume <- kernel.build(comm = t(presence), trait = traits)))
    
    ## Get the values (remove messages)
    suppressMessages(silent <- capture.output(richness <- BAT::kernel.alpha(comm=hypervolume )))
    suppressMessages(silent <- capture.output(dispersion <- BAT::kernel.dispersion(comm = hypervolume)))
    suppressMessages(silent <- capture.output(regularity <- BAT::kernel.evenness(comm = hypervolume)))

    ## Get the values
    return(cbind("hypervolume_richness"   = richness,
          "hypervolume_dispersion" = dispersion,
          "hypervolume_regularity" = regularity))
}

## Calculating the convex hull metrics
fun.convex.hull <- function(presence, traits, verbose) {
    if(verbose) {
        message(".", appendLF = FALSE)
    }

    ## get the hull
    suppressMessages(silent <- capture.output(hull <- hull.build(comm = t(presence), trait = traits)))
    
    # run the volume
    suppressMessages(silent <- capture.output(results <- BAT::hull.alpha(hull)))
    results <- t(t(results))
    rownames(results) <- colnames(presence)
    colnames(results) <- "convex.hull"
    return(results)
}

## Calculating the hypervolume metrics
fun.proba.den <- function(presence, traits, verbose) {
    if(verbose) {
        message(".", appendLF = FALSE)
    }

    ## Get the sds matrix
    standard_deviations <- matrix(
        data = rep(sqrt(diag(ks::Hpi.diag(traits))), nrow(presence)),
        nrow = nrow(presence),
        byrow = TRUE,
        dimnames = list(rownames(presence), c("sdA1", "sdA2")))

    ## Measure probability density metrics 
    suppressMessages(TPD_species <- TPD::TPDsMean(species = rownames(presence),
                                 means = traits,
                                 sds = standard_deviations,
                                 alpha = 0.99))
    TPD_cases   <- TPD::TPDc(TPDs = TPD_species, sampUnit = t(presence))
    output <- suppressMessages(do.call(cbind, TPD::REND(TPDc = TPD_cases)$communities))
    colnames(output) <- paste("TPD", gsub("F", "", colnames(output)), sep = "_")
    return(output)
}
