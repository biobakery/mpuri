#' Performs Mann-Witney U test on each group vs. control, each bug
#' 
#' \code{hypothesis_test} performs Mann-Witney U test, using
#' \code{\link{wilcox.test}} for each OTU in the data frame \code{d}. Each level
#' of \code{group}, excluding the \code{base_group} level, is tested against the
#' \code{base_group} level at an alpha level of 0.05 for evidence against the
#' null hypothesis that distributions of the the normalized counts for
#' \code{base_group} and the tested group do not differ by a location shift. The
#' alternative is two-sided.
#'
#' The data frame d is assumed to be in \emph{long} format. 
#'
#' @param d Data frame of data. 
#' @param group Name of data frame column containing the group information.
#' @param otuname_col Name of data frame column containing OTU names.
#' @param value_col Name of column containing normalized proportions for each of
#' the otunames.
#' @param base_group Name of the level in \code{group} indicating baseline
#' measurements.
#' @param alpha alpha level to conduct the Mann-Witney U test.
#' @return the same input data frame d, but with a new column called "changed".
#' This column is a factor with three levels: "no change" if the hypothesis test
#' was not significant, and either "decreased" or "increased" if the test was
#' significant. The sign of the change is calculated by the difference in means
#' for data from the tested group and data from the \code{base_group}.
hypothesis_test <- function(d, group, otuname_col = "bug", value_col = "value",
                            base_group = "Control", alpha = 0.05) {
    d$changed <- "no change"
    group_levels <- levels(d[[group]])
    for (bug in unique(d$bug)) {
        cur_bug <- d[d[[otuname_col]] == bug,]
        for (g in group_levels) { 
            if (g != base_group) {
                g1 <- cur_bug[cur_bug[[group]] == g, value_col]
                ctrl <- cur_bug[cur_bug[[group]] == base_group, value_col]
                res <- suppressWarnings( wilcox.test(g1, ctrl,
                                                     conf.level=(1-alpha)) )
                #res <- wilcox.test(g1, ctrl, conf.level=(1-alpha))
                if (is.na(res$p.value)) {
                    print("p-value is NA. Data is probably the same for both groups")
                    print(bug)
                    print(g)
                    print(g1)
                    print(ctrl)
                    print(res)
                    next
                }
                if (res$p.value <= alpha) {
                    #difference <- median(g1) - median(ctrl)
                    # should be using median but sometimes median is the same
                    # for both groups
                    difference <- mean(g1) - mean(ctrl)
                    if (difference < 0) {
                        d[d[[otuname_col]] == bug & d[[group]] == g, "changed"] <- "decreased" 
                    }
                    else {
                        if (difference > 0) {  
                            d[d[[otuname_col]] == bug & d[[group]] == g, "changed"] <- "increased" 
                        }
                        else {
                            print("Something weird has happened")
                            print(bug)
                            print(g)
                            print(g1)
                            print(ctrl)
                            print(res)
                        }
                    }
                }
            } 
        } # group for loop
    } # bug for loop
    d$changed <- factor(d$changed, levels=c("decreased", "no change",
                                            "increased"))
    return(d)
}

