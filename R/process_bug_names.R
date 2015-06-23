#' Pretty print taxonomy name.
#'
#' \code{process_bug_name} cleans up full taxonomy names to display the most
#' relevant information
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
process_bug_name <- function(name, sep=";") {

    if (stringr::str_sub(name, 1, 3) != "k__") {
        # silva style OTU name
        return(process_silva_name(name, sep))
    }
    else {
        # greengenes style OTU name
        return(process_gg_name(name, sep))
    }
}



#' Pretty print taxonomy name (greengenes-like).
#'
#' \code{process_bug_name} cleans up full taxonomy names to display the most
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
process_gg_name <- function(name, sep="; ") {
    result = ""
    for (level in class_levels) {
        # regex to detect 4 alphabetical characters. 
        # the [\\.]? in the beginning is because some had one period in
        # before the name. Not sure if that name was mis-coded or not
        reg <- paste0(level, "__[\\.]?[A-Za-z]{4,}[_A-Za-z0-9\\- ]*(", sep,
                    ")?")
        #print(reg)
        # detect some character. Used to display subspecies or
        # classification levels that just have a number
        reg_min <- paste0(level, "__[A-Za-z0-9_\\- ]+(", sep, ")?")
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


#' Pretty print taxonomy name (silva-like)
#'
#' \code{process_silva_name} cleans up full taxonomy names to display the most
#' relevant information.
#'
#' @param name Unprocessed, full taxonomy name. 
#' @param sep String separating different levels of taxonomic names
#' @return Return value will be a string that includes classification levels
#' with consecutive alphabetical characters as a sub-level of the main level <4
process_silva_name <- function(name, sep=";") {
    res <- stringr::str_split(name, sep)

    ind <- length(res[[1]])
    cur_name <- res[[1]][ind]
    result <- ""
    un_regex <- "[Uu]n(classified|identified|cultured).*"
    fullname_regex <- "[A-Za-z]{4,}[_A-Za-z0-9\\- ]*"
    for (ind in length(res[[1]]):1) {
        cur_name <- res[[1]][ind]
        is_unidentified <- !is.na(stringr::str_extract(cur_name, un_regex))
        is_fullname <- !is.na(stringr::str_extract(cur_name, fullname_regex))
        if (!is_unidentified && is_fullname) {
            result <- paste0(cur_name, ";", result)
            return(result)
        }
        else {
            result <- paste0(cur_name, ";", result)
        }
    }
    return(result)
}
