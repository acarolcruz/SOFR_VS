#' Check convergence of ELBO
#'
#' @param elbo_c 
#' @param elbo_prev 
#' @param convergence_threshold 
#'
#' @returns
#' @export
#'
#' @examples
check_convergence <- function(elbo_c, elbo_prev, convergence_threshold) {
  if(is.null(elbo_prev) == TRUE) {
    return(FALSE)
  }
  else{
    dif <- elbo_c - elbo_prev
    if(abs(dif)  <= convergence_threshold) return(TRUE)
    else return(FALSE)
  }
} 