#' @title paired.test
#'
#' @description paired t-test for a single stressor
#'
#' @param one_stressor the results from one stressor
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

paired.tests <- function(one_stressor) {
    ## Get the results table
    results_table <- one_stressor$results_table
    ## Split the results into two lists (random and stressor)
    random_list <- unlist(apply(results_table[, seq(from = 1, to = ncol(results_table), by = 2)], 2, list), recursive = FALSE)
    stressor_list <- unlist(apply(results_table[, seq(from = 2, to = ncol(results_table), by = 2)], 2, list), recursive = FALSE)
    ## Calculate all the differences
    test_list <- mapply(t.test, random_list, stressor_list, MoreArgs = list(paired = TRUE), SIMPLIFY = FALSE)
    ## Extract the parameters (d, t and p)
    extract.params <- function(one_test) {
        return(c(d = one_test$estimate[[1]],
                 t = one_test$statistic[[1]],
                 p = one_test$p.value[[1]]))
    }
    test_table <- t(do.call(cbind, lapply(test_list, extract.params)))
    rownames(test_table) <- gsub("random_", "", rownames(test_table))
    return(test_table)
}