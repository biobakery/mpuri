#' Plot boxplots of log10 relative abundance, facetted by 1 feature
#' 
#' \code{plot_group_panel} plots boxplots of log10 relative abundances, sorted
#' by their median and facetted by 1 feature. In each plot, the top
#' \code{cutoff} OTUs are plotted, as well as a line connecting the medians from
#' the \code{base_group}. Any increases or decreases in the other groups,
#' relative to the base group, are colored. The main inputs are a
#' \code{metadata} \code{\link{data.frame}}, and a matrix for the OTU table. The
#' metadata should be in \emph{long} format, with one row corresponding to one
#' sample. The OTU table should have the OTUs as the rows, and the samples as
#' the columns.
#'
#' \code{plot_group_panel} plots boxplots of log10 relative abundances, sorted
#' by their median and facetted by 1 feature. In each plot, the top
#' \code{cutoff} OTUs are plotted, as well as a line connecting the medians from
#' the \code{base_group}. Any increases or decreases in the other groups,
#' relative to the base group, are colored. The main inputs are a
#' \code{metadata} \code{\link{data.frame}}, and a matrix for the OTU table. The
#' metadata should be in \emph{long} format, with one row corresponding to one
#' sample. The OTU table should have the OTUs as the rows, and the samples as
#' the columns.
#'
#' The main inputs are a \code{metadata} \code{\link{data.frame}}, and a matrix
#' for the OTU table. The metadata should be in \emph{long} format, with one row
#' corresponding to one sample. The OTU table should have the OTUs as the rows,
#' and the samples as the columns.
#' 
#' The coloring is determined by the results of \code{\link{hypothesis_test}},
#' which currently uses a Mann-Witney U test as implemented in
#' \code{\link{wilcox.test}}.
#'
#' @param metadata Data frame of metadata, like the one returned by readPCL.
#' @param otu_tab_matrix matrix of relative abundance, like the one returned by
#' readPCL. Each row should be an OTU, and each column should be a sample.
#' @param group The grouping variable in the metadata.
#' @param base_group The reference level for the grouping variable
#' \code{group}. When performing hypothesis tests, all other levels of the
#' grouping variable \code{group} are compared to this one. 
#' @param cutoff Only display the top \code{cutoff} bugs in terms of median
#' relative abundance.
#' @param alpha alpha level to conduct each hypothesis test. 
#' @param notch Whether to have notched boxplots.
#' @param sep Character(s) separating different levels of taxonomic
#' classification in the OTU labels (rows of \code{otu_tab_matrix})
#' @param title Title for the graph.
#' @return ggplot2 plotting object. This function does not automatically display
#' the plot, so you need to call \code{print} on its output to display
#' the plot. 
#' @export
plot_group_panel <- function(metadata, otu_tab_matrix, group, 
                             base_group="Control",
                             cutoff=min(20, nrow(otu_tab_matrix)), 
                             alpha=0.05, notch=FALSE, sep=";", title="") {
    cleaned_data <- clean_plot_data(metadata, otu_tab_matrix, group, base_group,
                                    cutoff, alpha, sep)
    top_data <- cleaned_data[[1]]
    meds_plot <- cleaned_data[[2]]
    minval <- cleaned_data[[3]]
    maxval <- cleaned_data[[4]]
    meds_plot$type <- paste("Connecting medians of\ntaxa from", base_group)

    # change colors here for coloring no change/increase/decrease
    change_cols <- c("#0072B2", "#000000", "#D55E00")
    #change_cols <- c("#2C7BB6", "#FFFFBF", "#D7191C")
    names(change_cols) <- c("decreased", "no change", "increased")
    color_scale <- ggplot2::scale_colour_manual(name="change_cols",
                                                values=change_cols)
    # plotting code

    # update boxplot points so they can be the same color as the boxplot
    ggplot2::update_geom_defaults("point", list(colour=NULL))

    # draw some invisible points for the medians
    p1 <- ggplot2::ggplot(data=meds_plot, ggplot2::aes_string(x="bug_ord",
                                                              y="med")) +
            ggplot2::geom_point(alpha=0)

    # draw the boxplots
    p1 <- p1 + ggplot2::geom_boxplot(data=top_data, 
                                     ggplot2::aes_string(x="bug", y="value",
                                                         color="changed"), 
                                     notch=notch) 
    p1 <- p1 + color_scale
    p1 <- p1 + ggplot2::facet_grid(as.formula(paste("~", group)))
    # To change dots and line colors here:
    # draw the medians as points and connect them 
    p1 <- p1 + ggplot2::geom_point(ggplot2::aes_string(shape="type"),
                                   color="#999999") 
    p1 <- p1 + ggplot2::geom_line(ggplot2::aes_string(group="1",
                                                      linetype="type"),
                                  color="#999999")
    # To change X and Y label
    p1 <- p1 + ggplot2::xlab("Taxa") + ggplot2::ylab("log10 Relative Abundance")
    p1 <- p1 + ggplot2::coord_flip(ylim=c(minval, maxval))
    p1 <- p1 + ggplot2::ggtitle(title)
    p1 <- p1 + ggplot2::guides(colour = ggplot2::guide_legend(title = 
                                    sprintf("Significance: p<= %.4f", alpha)),
                               linetype = ggplot2::guide_legend(title=""),
                               shape = ggplot2::guide_legend(title=""))
    ggplot2::update_geom_defaults("point", list(colour="black"))
    return(p1)
}


clean_plot_data <- function(metadata, otu_tab_matrix, group, 
                            base_group="Control",
                            cutoff=min(20, nrow(otu_tab_matrix)), 
                            alpha=0.05, sep=";") {
    # get rid of R CMD check NOTE
    value <- NULL
    # get medians of control
    otu_tab.wide <- data.frame(t(otu_tab_matrix))
    otu_tab.wide[[group]] <- metadata[[group]]
    otu_tab.all <- reshape2::melt(otu_tab.wide, id=group, variable.name="bug")
    medians <- plyr::ddply(otu_tab.all, c("bug", group), plyr::summarise,
                           med=median(value))
    meds <- medians[medians[[group]]==base_group,]
    
    # sort medians and take the top cutoff
    orders <- order(meds$med)
    l <- length(orders)
    ord <- otu_tab_matrix[orders[(l-cutoff+1):l],]
    otu_tab.wide2 <- data.frame(t(ord), check.names=FALSE)
    
    # clean up names of bugs
    newnames <- sapply(colnames(otu_tab.wide2), 
                       function(x) { process_bug_name(x, sep=sep) })
    colnames(otu_tab.wide2) <- newnames
    # add grouping metadata
    otu_tab.wide2[[group]] <- metadata[[group]]
    
    # convert to long format
    d <- reshape2::melt(otu_tab.wide2, id=group, variable.name="bug") 
    
    # set truncation to deal with log10(0)
    minval <- log10(1e-16 + nzapply(d[, "value"], min)) 
    maxval <- log10(1e-16 + nzapply(d[, "value"], max))
    
    # log transform
    d$value <- log10(1e-16 + d$value)

    # test if significant
    d <- hypothesis_test(d, group, base_group=base_group, alpha=alpha)

    # get the medians
    meds_plot <- plyr::ddply(d, c("bug", group), plyr::summarise,
                             med=median(value), iqr=IQR(value), n=length(value))
    meds_plot <- meds_plot[meds_plot[[group]]==base_group,][, 
                                                c("bug", "med", "iqr", "n")]
    meds_plot$bug_ord <- with(meds_plot, 
                              reorder(bug, med, function(x) { x }, order=T))
    return(list(d, meds_plot, minval, maxval))
}


#' Plot everything in one panel
#'
#' \code{plot_one_panel} plots boxplots of log10 relative abundances, sorted by
#' their median. The top \code{cutoff} OTUs (in terms of relative abundance) are
#' plotted. The main inputs are a \code{metadata} \code{\link{data.frame}}, and
#' a matrix for the OTU table. The metadata should be in \emph{long} format,
#' with one row corresponding to one sample. The OTU table should have the OTUs
#' as the rows, and the samples as the columns.
#'
#' The main inputs are a \code{metadata} \code{\link{data.frame}}, and a matrix
#' for the OTU table. The metadata should be in \emph{long} format, with one row
#' corresponding to one sample. The OTU table should have the OTUs as the rows,
#' and the samples as the columns.
#' 
#' The coloring is determined by the results of \code{\link{hypothesis_test}},
#' which currently uses a Mann-Witney U test as implemented in
#' \code{\link{wilcox.test}}. 
#'
#' @param dots If \code{TRUE}, dots are plotted on top of the boxplots
#' @inheritParams plot_group_panel
#' @return ggplot2 plotting object. Must call print to actually display the
#' plot. 
#' @export
plot_one_panel <- function(metadata, otu_tab_matrix, group, 
                           base_group="Control", dots=TRUE,
                           cutoff=min(20, nrow(otu_tab_matrix)),
                           alpha=0.05, notch=FALSE, sep=";", title="") {
    cleaned_data <- clean_plot_data(metadata, otu_tab_matrix, group, base_group,
                                    cutoff, alpha, sep)
    top_data <- cleaned_data[[1]]
    meds_plot <- cleaned_data[[2]]
    minval <- cleaned_data[[3]]
    maxval <- cleaned_data[[4]]

    # plotting code
    # draw some invisible points for the medians
    p1 <- ggplot2::ggplot(data=meds_plot, ggplot2::aes_string(x="bug_ord",
                                                              y="med")) +
            ggplot2::geom_point(alpha=0)

    # draw the boxplots
    p1 <- p1 + ggplot2::geom_boxplot(data=top_data, 
                                     ggplot2::aes_string(x="bug", y="value",
                                                         color=group),
                                     notch=notch, outlier.shape=NA, 
                                     position="dodge")
    # draw rectangles to highlight certain regions
    # significant <- top_data[top_data$changed != "no change",]
    # rectangles <- data.frame(xmin=significant$bug, xmax=significant$bug,
    #                          ymin=-Inf, ymax=Inf)
    # p1 <- p1 + ggplot2::geom_rect(data=rectangles,
    #                               ggplot2::aes_string(xmin="xmin", xmax="xmax",
    #                                                   ymin="ymin", ymax="ymax"),
    #                               color="grey20", alpha=0.5, inherit.aes=FALSE)
    if (dots) {
        p1 <- p1 + ggplot2::geom_point(data=top_data,
                                       ggplot2::aes_string(x="bug", y="value",
                                                           color=group,
                                                           fill=group),
                                       position=ggplot2::position_jitterdodge(),
                                       alpha=0.5)
    }
    p1 <- p1 + ggplot2::xlab("Taxa") + ggplot2::ylab("log10 Relative Abundance")
    p1 <- p1 + ggplot2::coord_flip(ylim=c(minval, maxval))
    p1 <- p1 + ggplot2::ggtitle(title)
    return(p1)
}
