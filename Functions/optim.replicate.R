#' @title Optimised replicates
#'
#' @description A wrapper for replicating a based on results variance
#' 
#' @param fun The function to execute. The function should output a single value or a vector or a list of values.
#' @param diagnose The diagnosis function should intake a vector of values.
#' @param summarise Optional, if the output of \code{fun()} is complex, a function to summarise that output and pass it to \code{diagnose()}.
#' @param minimum The minimum of replicates to run (default is 20)
#' @param maximum Optional, a maximum of iterations (default is Inf).
#' @param stop.variance The percentage of variance change for stopping (default is 0.05).
#' @param increment The percentage of increment (if left empty, the value from stop.variance is used).
#' @param verbose Whether to be verbose (\code{TRUE}) or not (\code{FALSE}; default)
#' @param parallel Whether to use the \code{\link[future.apply]{future_replicate}} function for parallelisation (\code{TRUE}) or not (\code{FALSE}; default).
#' @param par.plan If \code{parallel = TRUE}, which plan to use (see \code{\link[future]{plan}}; by default, the plan is "future::multisession").
#' @param bkp.path Optional, a \code{"character"} string to be the path to a backup file. See details.
#' @param bkp.name Optional, if \code{bkp.path} is not missing, a \code{"character"} for the name of the file to save (or to load). See details.
#' @param ... Any optional arguments to be passed to the function.
#' 
#' @details When \code{bkp} is specified, the algorithm first looks if a file with \code{bkp.name} exists (by default, the name is \code{"optim.replicate.bkp"}). If not, the results will be saved after each series of replicates into this file. If it exist, the initial steps are skipped and the algorithm continues from the last backup point. \emph{NOTE} that this option slows down the simulations.
#' 
#' @return
#' A list of:
#' \itemize{
#'      \item{results}: a \code{matrix} with the summarised output from \code{fun()} for every set of iterations;
#'      \item{diagnoses}: a \code{matrix} with changes in variances after each diagnosis;
#'      \item{outputs}: a \code{list} of the outputs of each iteration of \code{fun()}.
#' }
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
optim.replicate <- function(fun, diagnose, summarise, minimum = 20, maximum = Inf, stop.variance = 0.05, increment, verbose = FALSE, parallel = FALSE, par.plan = "future::multisession", bkp.path, bkp.name, ...) {

    ## Check increment
    if(missing(increment)) {
        increment <- stop.variance
    }

    ## Check diagnosis
    if(missing(diagnose)) {
        stop("Impossible to optimize the replicate loop without a diagnose function.")
    }

    ## Check summary
    if(missing(summarise)) {
        summarise <- function(x) return(x)
    }

    ## Set the parallel (or not)
    if(parallel) {
        ## Set the plan for the future (please switch to: capitalism = NULL for a brighter future)
        future::plan(par.plan)
        ## Set the replicate function type
        rep.fun <- future.apply::future_replicate
    } else {
        rep.fun <- base::replicate
    }

    ## Set the bkp
    bkp_save <- bkp_exists <- FALSE
    if(!missing(bkp.path)) {
        ## Do the backup
        bkp_save <- TRUE
        if(missing(bkp.name)) {
            bkp.name <- "optim.replicate.bkp"
        }
        ## Does the backup file already exists?
        bkp_exists <- length(list.files(path = bkp.path, pattern = bkp.name)) > 0
    }

    if(bkp_exists) {
        ## Load the previous backup
        load(file = paste0(bkp.path, bkp.name, collapse = ""))
        n_iterations  <- bkp$n_iterations
        output_save   <- bkp$output_save
        results_table <- bkp$results_table
        diagnosis     <- bkp$diagnosis
        count <- n_iterations
    } else {
        ## Initialise the number of iterations
        count <- 0
    }

    ## Set the function to be verbose or not
    if(verbose) {
        if(!parallel) {
            ## Set the function to be verbose
            run.fun <- eval.verbose.count(fun, count = count)
        } else {
            run.fun <- eval.verbose(fun)
        }
    } else {
        run.fun <- fun
    }

    if(!bkp_exists) {
        ## Running the minimum number of iterations
        if(verbose) {
            message(paste0("Running the initial ", minimum, " replicates:" ), appendLF = FALSE)
        }
        output <- rep.fun(minimum, run.fun(), ..., simplify = FALSE)
        # output <- rep.fun(minimum, run.fun(), simplify = FALSE) ; warning("DEBUG")
        
        ## Saving the output and the number of iterations
        output_save <- output
        n_iterations <- minimum

        ## Summarizing the results
        results_table <- do.call(rbind, lapply(output, summarise))

        ## Iterations done
        if(verbose) {
            message("Done.", appendLF = FALSE)
        }

        ## Diagnosing the output variance
        diagnosis <- apply(results_table, 2, diagnose)
        diagnosis <- rbind(diagnosis, rep(0, length(diagnosis)))
        rownames(diagnosis) <- c(minimum, 0)

        if(bkp_save) {
            ## Save the progress
            bkp <- list("results_table" = results_table, "diagnosis" = diagnosis, "output_save" = output_save, "n_iterations" = n_iterations)
            save(bkp, file = paste0(bkp.path, bkp.name, collapse = ""))
        }
    }

    ## Continue the runs if any of the diagnoses has not reached the correct variance
    while(any(!abs(diagnosis[1,]/diagnosis[2,] - 1) < stop.variance) || n_iterations >= maximum) {
        
        ## Add some iterations
        add_iterations <- ceiling(n_iterations*increment)

        if(verbose) {
            message(paste0("\nRunning an additional ", add_iterations, " replicates:" ), appendLF = FALSE)
        }
        output <- rep.fun(add_iterations, run.fun(), ..., simplify = FALSE)
        # output <- rep.fun(add_iterations, run.fun(), simplify = FALSE) ; warning("DEBUG")

        ## Saving the output
        output_save <- c(output_save, output)
        n_iterations <- n_iterations + add_iterations

        ## Summarizing the results
        results_table <- rbind(results_table, do.call(rbind, lapply(output, summarise)))
        
        if(verbose) {
            message("Done.", appendLF = FALSE)
        }
    
        ## Diagnosing the output variance
        new_diagnosis <- apply(results_table, 2, diagnose)
        diagnosis <- rbind(new_diagnosis, diagnosis)
        rownames(diagnosis)[1] <- n_iterations

        if(verbose) {
            if(ncol(diagnosis) < 5) {
                message(paste0("\nDiagnosis change:" , paste0(round((diagnosis[2,]/diagnosis[1,]) - 1, digits = 3), collapse = ", ")), appendLF = FALSE)
            } else {
                ## Overall diagnosis change
                change <- (diagnosis[2,]/diagnosis[1,]) - 1
                round(mean(change), digits = 3)
                message(paste0("\nGlobal diagnosis change (n=", ncol(diagnosis), "): mean=", round(mean(change), digits = 3), "; sd=", round(sd(change), digits = 3), "."), appendLF = FALSE)
            }
        }

        if(bkp_save) {
            ## Save the progress
            bkp <- list("results_table" = results_table, "diagnosis" = diagnosis, "output_save" = output_save, "n_iterations" = n_iterations)
            save(bkp, file = paste0(bkp.path, bkp.name, collapse = ""))
        }        
    }

    if(verbose) {
        if(ncol(diagnosis) < 5) {
            message(paste0("\nResults converged after ", n_iterations, " iterations: additional variances (", paste0(round((diagnosis[2,]/diagnosis[1,]) - 1, digits = 3), collapse = ", "), ") all < ", stop.variance, "."), appendLF = FALSE)
        } else {
            ## Overall diagnosis change
            change <- (diagnosis[2,]/diagnosis[1,]) - 1
            round(mean(change), digits = 3)
            message(paste0("\nResults converged after ", n_iterations, " iterations: additional global variances (n=", ncol(diagnosis), "): mean=", round(mean(change), digits = 3), "; sd=", round(sd(change), digits = 3), ". All < ", stop.variance, "."), appendLF = FALSE)
        }
    }

    ## Return the output
    return(list("results" = results_table, "diagnoses" = diagnosis[rev(1:nrow(diagnosis)), , drop = FALSE], "outputs" = output_save))
}

eval.verbose <- function(fun, msg = ".") {
    new.fun <- function() {
        message(msg, appendLF = FALSE)
        return(fun())
    }
    return(new.fun)
}

eval.verbose.count <- function(fun, count) {
    new.fun <- function() {
        count <<- count + 1
        message(paste0(count, "."), appendLF = FALSE)
        return(fun())
    }
    return(new.fun)
}
