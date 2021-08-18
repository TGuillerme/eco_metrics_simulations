#' @title extract.table
#'
#' @description extract the scaled results per metric
#'
#' @param results the list of results (per stressor)
#' @param centre whether to centre the data on the null results (TRUE; default)
#' @param scale.method optional, which scaling method to use. Either "within" for scaling each metric within each stressor (abs(max(metric$stressor))), "between" for scaling each metric between each stressor (abs(max(unlist(metric)))) or left empty for no scaling.
#' @param reorder optional, change the order of the metrics (rows)
#' @param rename optional, change the names of the metrics
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

## Extract table
extract.table <- function(results, centre = TRUE, scale.method, reorder, rename) {


    # stop("DEBUG extract.table")
    # results <- all_results
    
    ## Centre all the results
    results <- lapply(results, centre.null, centre = centre)

    ## Scale all the metrics between each stressor
    if(!missing(scale.method)) {
        if(scale.method == "between") {
            ## Merging all the results together to the the max per metric
            metric_max <- apply(do.call(rbind, results), 2, function(x) return(max(abs(x))))
            ## Scaling each result by its maximum
            results <- lapply(results, function(X, max) t(apply(X, 1, function(x, max) return(x/max), max = metric_max)), max = metric_max)
        }
        ## Scale all the metrics within each stressor
        if(scale.method == "within") {
            ## Merging all the results together to the the max per metric
             results <- lapply(results, function(one_results) apply(one_results, 2, function(X) return(X/max(abs(X)))))
        }
    }

    ## Sort the results per metric
    results <- lapply(results, sort.results.names, centre = centre, reorder, rename)

    return(results)
}

## Centring on the null
centre.null <- function(one_results, centre = TRUE) {

    ## Finding the null results
    null_results <- grep("random", colnames(one_results$results_table))

    ## Get the pairs of columns
    metrics_names <- gsub("stressor_", "", gsub("random_", "", colnames(one_results$results_table)))
    pairs <- lapply(as.list(unique(metrics_names)), function(x, names) which(names %in% x), names = metrics_names)
    
    ## Centre each pair on the null
    if(centre) {
        centre_results <- lapply(pairs, function(pair, matrix) return(matrix[, pair[2]] - matrix[, pair[1]]), matrix = one_results$results_table)
        output <- do.call(cbind, centre_results)
    } else {
        output <- one_results$results_table
    }

    if(centre) {
        colnames(output) <- unique(metrics_names)
    } 
    return(output)
}

## Sort the results per name
sort.results.names <- function(results, centre, reorder, rename) {

    ## Find the different removal levels
    levels <- unique(gsub(".*_rm", "", colnames(results)))
     
    if(centre) {   
        ## Find the metrics names    
        metrics <- unique(gsub("_rm.*", "", colnames(results)))
    } else {
        ## Find the metrics names    
        metrics <- unique(gsub("stressor_", "", gsub("random_", "", gsub("_rm.*", "", colnames(results)))))     
    }

    if(!missing(reorder)) {
        metrics <- metrics[reorder]
    }

    ## Sort the metrics names
    sorted_metrics <- unlist(lapply(as.list(metrics), function(x, levels) paste(x, levels, sep = "_rm"), levels = levels))

    if(!centre) {
        ## Sort the metric names
        sorted_metrics <- c(apply(cbind(paste0("random_", sorted_metrics), paste0("stressor_", sorted_metrics)), 1, c))
    }

    ## Sort the column names
    results <- results[, match(sorted_metrics, colnames(results))]
    if(!missing(rename)) {
        ## TODO: optimise this bit (but just for the challenge)
        original_names <- metrics
        new_names <- rename
        col_names <- colnames(results)
        while(length(original_names) > 0) {
            col_names <- gsub(original_names[1], new_names[1], col_names)
            original_names <- original_names[-1]
            new_names <- new_names[-1]
        }
        colnames(results) <- col_names
    }
    ## Rename the column names
    return(results)
}

