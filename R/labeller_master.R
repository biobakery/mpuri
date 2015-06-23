#' Re-label facets when using ggplot2's \code{\link[ggplot2]{facet_grid}}.
#'
#' \code{labeller_master} returns a function that can be used to re-label factes
#' when using ggplot2's \code{\link[ggplot2]{facet_grid}}. 
#'
#' \code{labeller_master} returns a function that can be used to re-label factes
#' when using ggplot2's \code{\link[ggplot2]{facet_grid}}. The input list should
#' map the default facet labels to new facet labels. This function will only
#' work when facetting on 1 variable.
#
#' @param label_list A list mapping old labels to new labels.
#' @return A ggplot2 labeller function that can be passed as the \code{labeller}
#' argument for \code{\link[ggplot2]{facet_grid}}.
#'
#' @export
labeller_master <- function(label_list) {
    function(variable, value) {
        return(label_list[value])
    }
}
