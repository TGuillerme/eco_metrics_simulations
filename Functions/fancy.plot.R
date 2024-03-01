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
#' @param metric.names optional, a vector of fancy names for the results (rows)
#' @param ... plot options
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
fancy.plot <- function(results, probs = c(0.025, 0.25, 0.75, 0.975), cent.tend = median, col.metrics, lm, null, metric.names = NULL, ...) {

    ## Get the plot parameters
    plot_params <- list(...)
    if(is.null(plot_params$pch)) {
        plot_params$pch = 19
    }
    if(is.null(plot_params$lwd)) {
        plot_params$lwd = 1
    }

    if(length(results) > 1) {
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
                           grid  = FALSE,
                           metric.names = metric.names)

            } else {
                par(mar = c(5, 0, 4, 1))
                ## Making the empty plot
                empty.plot(results[[one_stressor]],
                           plot.main = names(results)[one_stressor],
                           add.x = TRUE,
                           add.y = FALSE)

                ## Adding the data
                if(!missing(lm)) {
                    model_data <- lm[[one_stressor]]
                } else {
                    model_data <- NULL
                }
                if(!missing(null)) {
                    null_data <- null[[one_stressor]]
                } else {
                    null_data <- NULL
                }
                ## Plot the data
                add.one.stressor(results[[one_stressor]], plot_params = plot_params, probs = probs, cent.tend = cent.tend, col = col.metrics, lm = model_data, null = null_data)
            }
        }
    } else {
        par(mfrow = c(1,1), mar = c(5, 10, 4, 1))
        ## Making the empty plot
        empty.plot(results[[1]],
                   plot.main = names(results)[1],
                   add.x = TRUE,
                   add.y = TRUE,
                   metric.names = metric.names)

        ## Adding the data
        if(!missing(lm)) {
            model_data <- lm[[1]]
        } else {
            model_data <- NULL
        }
        if(!missing(null)) {
            null_data <- null[[1]]
        } else {
            null_data <- NULL
        }
        ## Plot the data
        add.one.stressor(results[[1]], plot_params = plot_params, probs = probs, cent.tend = cent.tend, col = col.metrics, lm = model_data, null = null_data)
    }
}

## Adding one single line
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
    do.call(points, add_point)
}

## Adding the results for one metric
add.one.metric <- function(one_metric, plot_params, probs, cent.tend, col, lm_data = NULL, null_data = NULL) {
    
    ## Selecting the levels
    levels <- ncol(one_metric)

    ## Selecting the colours
    colours <- rev(sapply(seq(from = 0.2, to = 1, length.out = levels), function(x, col) adjustcolor(col, alpha.f = x), col = col))

    ## Get the plotting base
    base <- plot_params$y

    ## Add the linear model (if not null)
    if(!is.null(lm_data)) {
        ## Select the line significance
        lty <- ifelse(all(lm_data[c(2,4)] < 0.05), 1, 3)
        ## Set the coordinates of the lm
        y = c((1-0.2),(levels+0.2))
        ## Calculate the x coordinates (normally y = ax+b but flipped)
        x <- c(lm_data["Slope"]*y[2] + lm_data["Intercept"],
               lm_data["Slope"]*y[1] + lm_data["Intercept"])
        ## Adjust the y coordinates
        y = y+(base-levels)
        ## Add the model line
        lines(x, y, lty = lty, col = "darkgrey")

        ## Add the R^2
        text(paste(names(lm_data[5]), round(lm_data[5], 3), sep = "\n") , y = base-0.4, x = 0.6, cex = 1)
    }

    ## Plot each level
    for(one_level in 1:levels) {
        ## Setting the plot params
        plot_params$col <- colours[one_level]
        plot_params$y   <- rev(((base-levels)+1):base)[one_level]
        ## Setting the pch param
        if(!is.null(null_data)) {
            if(null_data[3*one_level] > 0.05) {
                plot_params$col <- "grey"
            }
            # plot_params$pch <- 1
            # if(null_data[3*one_level] < 0.05) {
            #     plot_params$pch <- 20
            #     if(null_data[3*one_level] < 0.001) {
            #         plot_params$pch <- 19
            #     }
            # }
        }
        add.one.line(one_metric[, one_level], plot_params = plot_params, probs = probs, cent.tend = cent.tend)
    }
}

## Adding the results for one stressor
add.one.stressor <- function(one_results, plot_params, probs, cent.tend, col, lm, null) {

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

        ## lm
        if(!is.null(lm)) {
            lm_data <- lm[i,]
        } else {
            lm_data <- NULL
        }

        ## null
        if(!is.null(null)) {
            null_data <- null[i,]
        } else {
            null_data <- NULL
        }

        ## Plotting the metric
        add.one.metric(one_results[, seq(from = 1 + (n_levels*(i-1)), to = i*n_levels, by = 1)], plot_params = plot_params, probs = probs, cent.tend = cent.tend, col = col[i], lm_data = lm_data, null_data = null_data)
    }
}

## Creating an empty plot
empty.plot <- function(one_results, scaled = TRUE, add.x = FALSE, add.y = FALSE, grid = TRUE,  plot.main, metric.names = NULL) {
    
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
    if(is.null(metric.names)) {
        metric_names <- fancy.names(split.metric(one_results, just.names = TRUE))
    } else {
        metric_names <- metric.names
    }
    n_levels <- ncol(one_results)/length(metric_names)


    ## Add the y labels
    if(add.y) {
        axis(2, at = 1:length(metric_names)*n_levels-median(1:n_levels)+1, labels = metric_names, las = 2)
    }

    if(grid) {
        ## Add the separators
        if(scaled) {
            ## Vertical 0 line
            abline(v = 0, col = "lightgrey")
        }
        ## Horizontal metric lines
        abline(h = 0:length(metric_names)*n_levels+0.5, col = "lightgrey", lty = 2)
    }
}