#' Pretty print taxonomy name
#'
#' \code{process_bug_names} cleans up full taxonomy names to display the most
#' relevant information
#' 
#' The classification levels in the name should be prefixed by k__, p__, etc.
#' for kingdom, phylum, etc (like greengenes profiles). 
#'
#' Tries to return genus + species, if species is present. If genus is not
#' present (e.g. contains less than 4 consecutive alphabetical characters),
#' tries to return the next highest level of classification. If that is not
#' present, keeps trying. 
#'
#' @param name Unprocessed, full taxonomy name. 
#' @param sep String separating different levels of taxonomic names
#' @return Return value will be a string that includes classification levels
#' with consecutive alphabetical characters as a sub-level of the main level <4
process_bug_names <- function(name, sep="; ") {
    result = ""
    for (level in class_levels) {
        # regex to detect 4 alphabetical characters. 
        # the [\\.]? in the beginning is because some had one period in before
        # the name. Not sure if that name was mis-coded or not
        reg <- paste0(level, "__[\\.]?[A-Za-z]{4,}[_A-Za-z0-9]*(", sep,
                      ")?")
        #print(reg)
        # detect some character. Used to display subspecies or classification
        # levels that just have a number
        reg_min <- paste0(level, "__[A-Za-z0-9_]+(", sep, ")?")
        m <- stringr::str_extract(name, reg)
        m_min <- stringr::str_extract(name, reg_min)
        if(!is.na(m)) {
            if (level == "s") {
                result <- m  
            }
            else {
               result <- paste0(m, result)
               return(result)
            }
        }
        else {
            if(!is.na(m_min)) {
                result <- paste0(m_min, result)
            }
        }
    }
    return(result)
}
