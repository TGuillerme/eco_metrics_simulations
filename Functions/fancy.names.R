#' @title Fancy names
#'
#' @description Makes fancy names out of columns of result tables (for plotting)
#'
#' @param names the column names
#' @param sep how to separate the underscores in names. Default is \n
#' @param remove.levels whether to remove eventual levels (e.g. _rm1). Default is TRUE.
#' @param nothing do nothing (just returns the non-fancied names). Default is FALSE.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

fancy.names <- function(names, sep = ":\n", remove.levels = TRUE, nothing = FALSE) {

    ## Do nothing
    if(nothing) {
        return(names)
    }

    ## Remove levels numbers
    if(remove.levels) {
        names <- gsub("_rm.", "", names)
    }
    ## Replace underscores
    return(gsub("_", sep, names))
}