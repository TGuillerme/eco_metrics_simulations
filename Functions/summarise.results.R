#' @title Summarise the results
#'
#' @description Summarises the results (from the replications)
#'
#' @param results the results or the path to the results
#' @param scale logical, whether to scale the results relative to the the random reduction (TRUE, default) or not (FALSE).
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

summarise.results <- function(results, scale = TRUE, cent.tend = median, prob = c(0.975, 0.75, 0.25, 0.025), output.by = "stressor") {

    ## Make the results into a list (or not)
    if(all(names(results) %in% c("results_table", "diagnosis", "output_save", "n_iterations"))) {
        results <- list(results)
    }

    ## Extracting the results table
    results_table <- extract.table(results, scale)

    ## Summarise the results
    return(lapply(results_table, summarise.table, cent.tend, prob))
}


## Centring on the null
centre.null.scale <- function(results, scale = FALSE) {

    ## Finding the null results
    null_results <- grep("random", colnames(results$results_table))

    ## Get the pairs of columns
    metrics_names <- gsub("stressor_", "", gsub("random_", "", colnames(results$results_table)))
    pairs <- lapply(as.list(unique(metrics_names)), function(x, names) which(names %in% x), names = metrics_names)
    
    ## Centre each pair on the null
    centre_results <- lapply(pairs, function(pair, matrix) return(matrix[, pair[2]] - matrix[, pair[1]]), matrix = results$results_table)
    output <- do.call(cbind, centre_results)

    ## Scale each centred results
    if(scale) {
        output <- apply(output, 2, function(X) return(X/max(abs(X))))
    }
    colnames(output) <- unique(metrics_names)
    return(output)
}

## Sort the results per name
sort.results.names <- function(results) {
    ## Find the different removal levels
    levels <- unique(gsub(".*_rm", "", colnames(results)))
    ## Find the metrics names    
    metrics <- unique(gsub("_rm.*", "", colnames(results)))
    ## Sort the metrics names
    sorted_metrics <- unlist(lapply(as.list(metrics), function(x, levels) paste(x, levels, sep = "_rm"), levels = levels))
    ## Return the sorted column names
    return(results[, match(sorted_metrics, colnames(results))])
}

## Extract table
extract.table <- function(results, scale) {
    ## Centre and scale the results
    results <- lapply(results, centre.null.scale, scale = scale)

    ## Sort the results per metric
    results <- lapply(results, sort.results.names)    

    return(results)
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