#' @title Summarise the results
#'
#' @description Summarises the results (from the replications)
#'
#' @param results the results or the path to the results
#' @param scale.method optional, which scaling method to use. Either "within" for scaling each metric within each stressor (abs(max(metric$stressor))), "between" for scaling each metric between each stressor (abs(max(unlist(metric)))) or left empty for no scaling.
#' @param centre whether to centre the data on the null results (TRUE; default)
#' @param cent.tend function, the central tendency function (default is median)
#' @param prob the probability quantiles to display (default is c(0.975, 0.75, 0.25, 0.025))
#' @param output.by whether to return the output sorted by stressor ("stressor"; default) or by metric ("metric").
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

summarise.results <- function(results, scale.method = "between", centre = TRUE, cent.tend = median, prob = c(0.975, 0.75, 0.25, 0.025), output.by = "stressor") {

    ## Make the results into a list (or not)
    if(all(names(results) %in% c("results_table", "diagnosis", "output_save", "n_iterations"))) {
        results <- list(results)
    }

    ## Extracting the results table
    results_table <- extract.table(results, centre, scale.method)

    ## Summarise the results
    return(lapply(results_table, summarise.table, cent.tend, prob))
}


## Summarise the results
summarise.table <- function(results, cent.tend, prob) {
    ## Central tendencies
    cent_tends <- apply(results, 2, cent.tend)
    ## Quantiles
    quants <- t(apply(results, 2, quantile, probs = prob))
    ## Get the metric names
    metric_names <- gsub("_rm.", "", colnames(results))
    ## Removal
    remove_level <- gsub(".*_rm", "", colnames(results))
    ## Combining everything
    table_out <- cbind("metric" = metric_names, "removal level" = remove_level)
    values_out <- cbind(cent_tends, quants)
    return(data.frame(cbind(table_out, values_out)))
}