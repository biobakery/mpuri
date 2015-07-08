#' Apply a vectorized function to all nonzero elements of a vector.
#'
#' \code{nzapply} applies a vectorized function to all nonzero elements of a
#' vector.
#'
#' Source:
#'
#'   \url{http://stackoverflow.com/questions/12583569}
#'
#'
#' @param x A vector.
#' @param FUN A vectorized function.
#' @param base What to return if all elements of x are zero.
#' @param ... Other arguments to pass to \code{FUN}
#' @return Returns the value of applying \code{FUN} to nonzero elements of
#'   \code{x}
#'
#' @examples
#' mpuri:::nzapply(c(0, 1, 2, 3, 4), min)
#' mpuri:::nzapply(c(0, 0, 0, 0), max, base = NA)
nzapply <- function(x, FUN, base=0, ...) {
    # apply a function func to the nonzero values of x
    # x: a vector
    # func: a function that takes a vector
    # base: base case when all the values of x are 0
    # ...: additional arguments to pass to func
    zerovals <- (x==0)
    if (all(zerovals)) {
        return(base)
    } else {
        return( FUN(x[!zerovals], ...) )
    }
}
