#' @title Plot results
#'
#' @description Plotting the results
#'
#' @param results the different results output from 02_Simulation_pipeline.Rmd
#' @param metric which metric to display?
#' @param scale.method optional, which scaling method to use. Either "within" for scaling each metric within each stressor (abs(max(metric$stressor))), "between" for scaling each metric between each stressor (abs(max(unlist(metric)))) or left empty for no scaling.
#' @param centre whether to centre all the metrics on the random (null) changes (TRUE; default). If set to FALSE the results are not centred and the null changes are also displayed.
#' @param level.names which level names to use (can be left empty)
#' @param legend whether to add the legend (TRUE, default) or not (FALSE)
#' @param legend.pos where to plot the legend (default is "topright")
#' @param ... any arguments to be passed to \code{boxplot}.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
plot.metric.results <- function(results, metric, scale.method = "between", centre = TRUE, level.names, legend = TRUE, legend.pos = "topright", ...) {

    ## Make the results into a list (or not)
    if(!is.null(names(results)) && all(names(results) %in% c("results_table", "diagnosis", "output_save", "n_iterations"))) {
        results <- list(results)
    }

    ## Extracting the results table
    results_table <- extract.table(results, centre, scale.method = scale.method)
    names(results_table) <- names(results)

    ## Plot the metric results
    plot.one.metric(results_table, metric, scale.method, centre, level.names, legend, legend.pos,...)
}

## Plotting one series of results
plot.one.metric <- function(results_table, metric, scale.method, centre, level.names, legend, legend.pos, ...) {

    plot_args <- list(...)

    ## Select the correct metric
    plot_data <- lapply(results_table, function(x, metric) x[, gsub("_rm.", "", colnames(x)) %in% metric], metric = metric)

    ## Trim lengths
    plot_data <- lapply(plot_data, function(x, min) x[sample(1:nrow(x), min), ], min = min(unlist(lapply(plot_data, nrow))))

    ## Combine plots
    plot_data <- do.call(cbind, plot_data)

    ## Get the number of levels
    levels <- unique(gsub(".*_rm", "", colnames(results_table[[1]])))

    ## Adding the level names
    if(missing(level.names)) {
        level.names <- levels
    }
    colnames(plot_data) <- rep(level.names, ncol(plot_data)/length(levels))
    
    plot_args$x <- plot_data

    ## Graphics handling
    ## y limit
    if(is.null(plot_args$ylim)) {
        if(!missing(scale.method)) {
            plot_args$ylim <- c(-1, 1)
        } else {
            plot_args$ylim <- range(plot_data)
        }
    }

    ## colour
    if(is.null(plot_args$col)) {
        plot_args$col <- rainbow(ncol(plot_data)/length(levels))
    } else {
        plot_args$col <- plot_args$col
    }
    ## Create the colour gradients
    legend_col <- plot_args$col
    plot_args$col <- unlist(lapply(as.list(plot_args$col), function(cols, levels) return(rev(sapply(seq(from = 0.2, to = 1, length.out = length(levels)), function(x, col) adjustcolor(col, alpha.f = x), col = cols))), levels))

    ## main
    if(is.null(plot_args$main)) {
        plot_args$main <- metric
    }

    ## Plot the stuff
    do.call(boxplot, plot_args)
    abline(h = 0, col = "lightgrey")

    if(legend) {
        legend(legend.pos, legend = names(results_table), col = legend_col, pch = 15)
    }
}
