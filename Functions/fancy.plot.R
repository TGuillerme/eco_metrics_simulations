#' @title Fancy plot
#'
#' @description Makes the fancy global plots
#'
#' @param results all results
#' @param probs the quantiles to display (default is c(0.025, 0.25, 0.75, 0.975))
#' @param cent.tend the central tendency function (default is median)
#' @param col.metrics the colours for each metric
#' @param lm optional, lm results to add to the plot
#' @param null optional, paired t-tests to add to the plot
#' @param ... plot options
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
fancy.plot <- function(results, probs = c(0.025, 0.25, 0.75, 0.975), cent.tend = median, col.metrics, lm, null, ...) {

    ## Get the plot parameters
    plot_params <- list(...)
    if(is.null(plot_params$pch)) {
        plot_params$pch = 19
    }
    if(is.null(plot_params$lwd)) {
        plot_params$lwd = 1
    }

    ## Setting up the main plotting window
    par(mfrow = c(1, length(results)+1))
    ## Plotting each stressor
    for(one_stressor in 0:length(results)) {
        
        if(one_stressor == 0) {
            ## Making the labels plot
            par(mar = c(5, 10, 4, 0), bty = "n")
            empty.plot(results[[1]],
                       plot.main = "",
                       add.x = FALSE,
                       add.y = TRUE,
                       grid  = FALSE)

        } else {
            par(mar = c(5, 0, 4, 1))
            ## Making the empty plot
            empty.plot(results[[one_stressor]],
                       plot.main = names(results)[one_stressor],
                       add.x = TRUE,
                       add.y = FALSE)

            ## Adding the data
            add.one.stressor(results[[one_stressor]], plot_params = plot_params, probs = probs, cent.tend = cent.tend, col = col.metrics)
        }
    }
}

## Adding one single line
#plot(NULL, xlim = c(-1, 1), ylim = c(0, 5))
#one_line <- one_metric[, 1]
#plot_params <- list(y = 1, lwd = 1, col = "blue", pch = 19)
#cent.tend = median
#probs = c(0.025, 0.25, 0.75, 0.975)
add.one.line <- function(one_line, plot_params, probs, cent.tend) {

    ## Get the number of quantiles
    n_quantiles <- (length(probs)/2)

    ## Measuring the quantiles
    ci_data <- quantile(one_line, probs = probs)

    ## Plotting each line
    for(cis in 1:n_quantiles) {
        ## Setting the line parameters
        add_line <- plot_params
        add_line$y <- rep(plot_params$y, 2)
        ## Getting the line width
        add_line$lwd <- plot_params$lwd + ((cis-1) * 1.5)
        ## Getting the line type
        add_line$lty <- n_quantiles - (cis - 1)
        ## Getting the line data
        add_line$x <- ci_data[c(cis, n_quantiles*2 - (cis-1))]
        ## Add the line
        do.call(lines, add_line)
    }

    ## Adding the median
    add_point <- plot_params
    add_point$x <- cent.tend(one_line)
    add_point$pch <- 19
    do.call(points, add_point)
}

## Adding the results for one metric
#one_metric <- one_results_sub[, 1:4]
add.one.metric <- function(one_metric, plot_params, probs, cent.tend, col) {
    
    ## Selecting the levels
    levels <- ncol(one_metric)

    ## Selecting the colours
    colours <- rev(sapply(seq(from = 0.2, to = 1, length.out = levels), function(x, col) adjustcolor(col, alpha.f = x), col = col))

    ## Get the plotting base
    base <- plot_params$y

    ## Plot each level
    for(one_level in 1:levels) {
        ## Setting the plot params
        plot_params$col <- colours[one_level]
        plot_params$y   <- rev(((base-levels)+1):base)[one_level]
        add.one.line(one_metric[, one_level], plot_params = plot_params, probs = probs, cent.tend = cent.tend)
    }
}

## Adding the results for one stressor
#one_results <- results[[1]]
add.one.stressor <- function(one_results, plot_params, probs, cent.tend, col) {

    ## Number of metrics
    n_metrics <- length(split.metric(one_results, just.names = TRUE))
    ## Number of levels
    n_levels <- ncol(one_results)/n_metrics

    ## Set the colours
    if(length(col) != n_metrics) {
        col <- rep(col, n_metrics)
    }

    ## Adding each metric
    for(i in 1:n_metrics) {
        ## Update the plotting location
        plot_params$y <- rev(1:n_metrics*n_levels)[i]
        ## Plotting the metric
        add.one.metric(one_results[, seq(from = 1 + (n_levels*(i-1)), to = i*n_levels, by = 1)], plot_params = plot_params, probs = probs, cent.tend = cent.tend, col = col[i])
    }
}

## Creating an empty plot
empty.plot <- function(one_results, scaled = TRUE, add.x = FALSE, add.y = FALSE, grid = TRUE, plot.main) {
    
    ## Number of metrics
    y_lim <- c(1,ncol(one_results))
    
    ## Scale values
    if(scaled) {
        x_lim <- c(-1, 1)
    } else {
        x_lim <- range(one_results)
    }
    
    ## Base plot
    if(add.x) {
        plot(NULL, xlim = x_lim, ylim = y_lim, ylab = "", xlab = "scaled metric change", yaxt = "n", main = plot.main)
    } else {
        plot(NULL, xlim = x_lim, ylim = y_lim, ylab = "", xlab = "", xaxt = "n", yaxt = "n", main = plot.main)
    }

    ## Get the metric names
    metric_names <- split.metric(one_results, just.names = TRUE)

    ## Add the y labels
    if(add.y) {
        n_levels <- ncol(one_results)/length(metric_names)
        axis(2, at = 1:length(metric_names)*n_levels-median(1:n_levels)+1, labels = metric_names, las = 2)
    }

    if(grid) {
        ## Add the separators
        if(scaled) {
            ## Vertical 0 line
            abline(v = 0, col = "grey")
        }
        ## Horizontal metric lines
        abline(h = 0:length(metric_names)*4+0.5, col = "grey", lty = 2)
    }
}