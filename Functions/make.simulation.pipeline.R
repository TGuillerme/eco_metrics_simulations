#' @title Make the simulation pipeline
#'
#' @description Set the variables for running one simulation
#'
#' @param type Which removal type (either \code{"facilitation"}, \code{"equalizing"}, \code{"competition"} or \code{"filtering"})}) to be paired with \code{"random"}.
#' @param sim.data A named list of parameters for the data simulation (default is \code{list{n.traits = 2, speciation = 1, n.taxa = 200}})
#' @param remove How much data to remove (can be a vector of multiple percentages)
#' @param verbose whether to be verbose (\code{TRUE}, default) or not (\code{FALSE})
#' @param record.timer whether to record the timer for each step of the simulations#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

make.simulation.pipeline <- function(type = NULL,
                                     sim.data = list(n.traits   = 2,
                                                     speciation = 1,
                                                     n.taxa     = 200),
                                     remove = c(0.2, 0.4, 0.6, 0.8),
                                     verbose  = TRUE,
                                     record.timer = TRUE) {
    return(function() simulation.pipeline(sim.data, type, remove, verbose = verbose, record.timer = record.timer))
}