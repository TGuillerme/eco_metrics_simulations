#' @title fit.model.rm
#'
#' @description Fit a model on the different levels of removal per metric
#'
#' @param one_stressor the results from one stressor
#' @param model.fun the type of model function (kind of only works for lm for now)
#' @param formula the model (default is metric ~ rm for the metric value function of the removal level)
#' @param summary.out whether to summarise the output (TRUE; default; returns the intercept, slope, adjust R2 and their p.values) or not (FALSE)
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
fit.model.rm <- function(one_stressor, model.fun = lm, formula = metric ~ rm, summary.out = TRUE) {
    ## Get the data per metric
    metric_data <- split.metric(one_stressor)

    ## Table the data
    table.data <- function(X) {
        n_rows <- nrow(X)
        return(data.frame("rm" = rep(1:ncol(X), each = n_rows), "metric" = c(X)))
    }
    metric_data <- lapply(metric_data, table.data)

    ## Apply the models
    all_models <- lapply(metric_data, function(data, fun, formula) fun(formula, data), fun = model.fun, formula = formula)
    
    if(!summary.out) {
        return(all_models)
    } else {
        return(t(do.call(cbind, lapply(all_models, summarise.models))))
    }

    # res <- t(do.call(cbind, lapply(all_models, summarise.models)))
    # for(i in 1:15) {
    #     plot(metric_data[[i]])
    #     abline(a = res[i,1], b = res[i,3])
    #     Sys.sleep(2)
    # }
}

## Extract the summary
summarise.models <- function(one_model){
    model <- summary(one_model)
    results <- c(model$coefficients[1,c(1,4)], model$coefficients[2,c(1,4)], "adj.R^2" = model$adj.r.squared)
    names(results)[c(1,3)] <- c("Intercept", "Slope")
    return(results)

}