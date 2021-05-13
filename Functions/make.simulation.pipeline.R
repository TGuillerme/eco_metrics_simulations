#' @title Make the simulation pipeline
#'
#' @description Set the variables for running one simulation
#'
#' @param sim.data A named list of parameters for the data simulation (default is \code{list{n.traits = 2, speciation = 1, n.taxa = 200}})
#' @param verbose whether to be verbose (\code{TRUE}, default) or not (\code{FALSE})
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

make.simulation.pipeline <- function(type,
                                     sim.data = list(n.traits   = 2,
                                                     speciation = 1,
                                                     n.taxa     = 200),
                                     remove = c(0.2, 0.4, 0.6, 0.8),
                                     verbose  = TRUE,
                                     record.timer = TRUE) {
    return(function() simulation.pipeline(sim.data, type, remove, verbose = verbose, record.timer = record.timer))
}



n_traits <- 2
speciation_rate <- 1
n_taxa <- 200

    return()
}