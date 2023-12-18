#' @title Testing a metric on a space
#'
#' @description Testing whether a metric can pick up changes in a traitspace following a stressor.
#TG: in other words, this basically redo each cells in the table/Figure 2 in the manuscript (each cell being one stressor + on metric)
#'
#' @param data either the full traitspace (a matrix provided by the user) or a list of arguments to d
#' @param stressor either one of the stressors from dispRity::reduce.space or a list of TRUE and FALSE corresponding to each colum
#TG: for example for the observed bird extinctions, the list can be something like:
#my_stressor <- list("extinction_18thC" = c(TRUE, TRUE, TRUE, FALSE, FALSE), # meaning the first three species are extinct and the fourth and the fifth are not (with the species number corresponding to their rowname position in data) 
#                    "extinction_19thC" = c(TRUE, TRUE, FALSE, FALSE, FALSE)) 
#' @param metric the metric to apply to the data and stressor
#' @param replicates (for simulations only)

test.metric.space <- function(data, stressor, metric, replicates) {
    return(NULL)
}
summary.test.metric.space <- function(x, ...) {
    return(NULL)
}
plot.test.metric.space <- function(x, ...) {
    return(NULL)
}

