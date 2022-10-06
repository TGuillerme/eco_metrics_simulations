## Function for simplifying the plots
#' @param data the dataset
#' @param reduction the vector of selected element for the reduction
#' @param null the vector of selected element for the null
#' @param col two colours
#' @param pch the pch (default = 19)
#' @param trait which trait to plot (default is 1 for plot.1D and c(1,2) for plot.2D)
#' @param xlab the xlab
#' @param legend the legend text
plot.1D <- function(data, reduction, null, main, col = c("blue", "orange"), pch = 19, trait = 1, xlab = "trait values", legend) {
    ## Get the histograms
    main_plot <- hist(data, plot = FALSE)
    plot_1 <- hist(data[null, trait], plot = FALSE, breaks = main_plot$breaks)
    plot_2 <- hist(data[reduction, trait], plot = FALSE, breaks = main_plot$breaks)
    plot_height <- max(plot_1$counts, plot_2$counts)

    ## Make the colours as transparent (for overlaps)
    make.transparent <- function(col, alpha = 0.5) {
      rgb_vals <- grDevices::col2rgb(col)
      rgb_args <- list(red   = rgb_vals["red", ],
                       green = rgb_vals["green", ],
                       blue  = rgb_vals["blue", ],
                       maxColorValue = 255,
                       alpha = alpha * 255)
      return(do.call(grDevices::rgb, rgb_args))
    }
    ## Plot the histograms
    plot(plot_1, col = make.transparent(col[1]), ylim = c(0, plot_height),
        main = ifelse(missing(main), "", main),
        xlab = "Trait 1", ylab = "Counts")
    plot(plot_2, col = make.transparent(col[2]), add = TRUE)

    ## Adding the legend
    if(!missing(legend)) {
        legend(x = "topleft", col = col, pch = pch, legend = legend)
    }
    ## Add the points
    points(x = data[null, trait], y = rep(0, sum(null)), col = "black", bg = col[1], pch = 21)
    points(x = data[reduction, trait], y = rep(-plot_height/50, sum(reduction)), col = "black", bg = col[2], pch = 21)
}
plot.2D <- function(data, reduction, null, main, col = c("blue", "orange"), pch= 19, trait = c(1,2), legend, ...) {
    ## Plotting the histogram
    plot(NULL, xlim = range(trait_space[, trait[1]]), ylim = range(trait_space[, trait[2]]), main = main, xlab = "Trait 1", ylab = "Trait 2")
    points(trait_space[null, trait], pch = pch, col = col[1])
    points(trait_space[reduction, trait], pch = pch, col = col[2])
    ## Adding the legend
    if(!missing(legend)) {
        legend(x = "topleft", col = col, pch = pch, legend = legend)
    }
}