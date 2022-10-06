#' @title Analyse replicates
#' 
#' @description Analyse the replicates results
#'
#' @param results the results list or a path to the file
#' @param what optional, a sub-selection of the results to plot (either the names of the indices)
#' @param scale whether to scale the variance for each diagnoses (TRUE) or not (FALSE, default)
#' @param net.change whether to plot the variance net change (TRUE, default) or not (FALSE).
#' @param ... any arguments to be passed to \code{lines}.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

analyse.replicates <- function(results, what, scale = FALSE, net.change = TRUE, ...) {

    ## Loading the results from a backup file
    if(is(results, "character")) {
        name <- load(results)
        eval(parse(text = paste0("results <- ", name)))
    }

    ## Extracting the diagnoses
    diagnoses <- results$diagnosis

    ## What to plot
    if(missing(what)) {
        what <- colnames(results$diagnosis)
    } else {
        if(is(what, "numeric")) {
            what <- colnames(results$diagnosis)[what]
        }
    }

    ## Getting the net changes
    if(net.change) {
        diagnoses <- abs(diagnoses[-nrow(diagnoses), ]/diagnoses[-1, ]-1)
        diagnoses[nrow(diagnoses), ] <- 0
    }

    ## Scaling the values
    if(scale) {
        diagnoses <- apply(diagnoses, 2, function(x) x/max(x))
    }

    ## plot_params
    plot_params <- list(...)
    if(is.null(plot_params$xlim)) {
        plot_params$xlim <- range(as.numeric(rownames(diagnoses)))
    }
    if(is.null(plot_params$ylim)) {
        plot_params$ylim <- range(diagnoses)
    }
    if(is.null(plot_params$ylab)) {
        plot_params$ylab <- "variance"
        if(net.change) {
            plot_params$ylab <- paste0("net ", plot_params$ylab, " change")
        }
        if(scale) {
            plot_params$ylab <- paste0("scaled ", plot_params$ylab)
        }
    }
    if(is.null(plot_params$xlab)) {
        plot_params$xlab <- "iterations"
    }
    if(is.null(plot_params$lty)) {
        plot_params$lty <- 1
    }
    if(is.null(plot_params$lwd)) {
        plot_params$lwd <- 1
    }

    ## Setting the colours
    if(is.null(plot_params$col)) {
        col_vector <- grDevices::rainbow(length(what))
    } else {
        col_vector <- plot_params$col
        plot_params$col <- NULL
    }

    ## Adding the empty plot
    empty_plot <- list(x = NULL)
    empty_plot <- c(empty_plot, plot_params)
    do.call(plot, empty_plot)

    ## Adding all the lines
    add.one.line <- function(data, plot_params, col) {
        line_param <- plot_params
        line_param$x <- as.numeric(names(data))
        line_param$y <- data
        line_param$col <- col
        do.call(lines, line_param)
    }
    for(one_line in 1:length(what)) {
        add.one.line(diagnoses[, what[one_line]], plot_params, col = col_vector[one_line])
    }
    return(invisible())
}
