#' @title The diagnosis function for optim.replicate
#'
#' @description 
#'
#' @param 
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

optim.diagnose <- function() {

    ## Combine all the replicates together
    #test
    return()
}

#' @title The summary function for optim.replicate
#'
#' @description 
#'
#' @param 
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
#' 
optim.summarise <- function(x) {
    ## Get all the values
    output <- c(unlist(x$results))

    ## Names for one result
    names <- c(paste(rownames(x$results[[1]])[1], colnames(x$results[[1]]), sep = "_"),
               paste(rownames(x$results[[1]])[2], colnames(x$results[[1]]), sep = "_"))
    names <- paste(names, rep(paste0("rm", rep(1:length(x$results))), each = length(names)), sep = "_")

    names(output) <- names

    return(output)
}