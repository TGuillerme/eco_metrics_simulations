#' @title fancy.table
#'
#' @description makes a fancy table per metric and stressor
#'
#' @param data the list of results per stressor from which to get the columns
#' @param columns the list of columns or pair of columns. For pairs of columns, the second column should be the p-value and is replaced by a symbol.
#' @param symbols a vector of up to 5 symbols per p-value levels (see details)
#' @param digits the number of digits for rounding (default is 3)
#' 
#' @details
#' By default the symbols are c(" ", ".", "*", "**", "***") for 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

# data[[1]][c(8,9),  c(1, 4)] #OK
# extracted_data[[1]][c(8,9), , drop = FALSE] #OK

fancy.table <- function(data, columns, symbols = c("", ".", "*", "**", "***"), digits = 2) {
    ## Making the fancy table
    return(do.call(cbind, lapply(data, get.columns, columns, symbols, digits)))
}

## Get the fancy columns for one stressor
get.columns <- function(one_stressor, columns, symbols, digits) {
    ## Extracting the columns
    extracted_data <- lapply(columns, function(col, dat) dat[, col, drop = FALSE], dat = one_stressor)
    ## Fancied data
    fancy_data <- lapply(extracted_data, function(col, symbols, digits) apply(col, 1, p.symbol, symbols, digits), symbols = symbols, digits = digits)

    ## Combine that
    fancy_data <- as.data.frame(do.call(cbind, fancy_data))
    colnames(fancy_data) <- colnames(one_stressor)[unlist(lapply(columns, function(x) x[1]))]
    return(fancy_data)
}

## Get the symbol for each value
p.symbol <- function(x, symbols, digits) {
    ## Arbitrary p_value limits
    limits <- c(0.1, 0.05, 0.01, 0.001)

    if(length(x) == 1) {
        ## No p-value symbol needed
        return(round(x, digits))
    } else {
        ## Adding the p symbol
        if(x[2] > limits[1]) {
            return(paste0(round(x[1], digits), symbols[1]))
        } else {
            if(x[2] > limits[2]) {
                return(paste0(round(x[1], digits), symbols[2]))
            } else {
                if(x[2] > limits[3]) {
                    return(paste0(round(x[1], digits), symbols[3]))
                } else {
                    if(x[2] > limits[4]) {
                        return(paste0(round(x[1], digits), symbols[4]))
                    } else {
                        if(x[2] <= limits[4]) {
                            return(paste0(round(x[1], digits), symbols[5]))
                        }
                    }
                }
            }
        }
    }
}