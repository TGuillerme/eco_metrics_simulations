#' @title split.metric
#'
#' @description splits a table of metrics (typically from extract.table of results$results_table)
#'
#' @param table the table of metrics results
#' @param stressor whether to split by stressor/random (TRUE) or not (FALSE; default)
#' @param just.names whether to return just the metric names (TRUE) or not (FALSE; default)
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

split.metric <- function(table, stressor = FALSE, just.names = FALSE) {

    ## Get all the metric names
    metric_names <- unique(gsub("_rm.", "", colnames(table)))

    ## Split by stressor (if required)
    if(stressor) {
        metric_names <- unique(gsub("stressor_", "", gsub("random_", "", metric_names)))
    }

    ## Return just the names?
    if(just.names) {
        return(metric_names)
    }

    ## Get the list of splits
    splits <- lapply(as.list(metric_names), function(X, table) grep(X, colnames(table)), table = table)

    ## Splitting the data
    return(lapply(splits, function(split, table) return(table[, split]), table = table))
}