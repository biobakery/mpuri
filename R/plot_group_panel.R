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
#' @param title Title for the graph.
#' @param labeller labeller function to pass to ggplot2's
#' \code{\link[ggplot2]{facet_grid}}. The default uses the values of the labels. 
#' @return ggplot2 plotting object. This function does not automatically display
#' the plot, so you need to call \code{\link{print}} on its output to display
#' the plot. 
plot_group_panel <- function(metadata, otu_tab_matrix, group, 
                             base_group="Control",
                             cutoff=min(20, nrow(otu_tab_matrix)), 
                             alpha=0.05,
                             notch=F, title="", labeller="label_value") {
    
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
    otu_tab.wide2 <- data.frame(t(ord))
    
    # clean up names of bugs
    newnames <- sapply(colnames(otu_tab.wide2), process_bug_names) 
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
    d <- hypothesis_test(d, group, base_group=base_group, alpha=alpha)

    # get the medians
    meds_plot <- plyr::ddply(d, c("bug", group), plyr::summarise,
                             med=median(value), iqr=IQR(value), n=length(value))
    meds_plot <- meds_plot[meds_plot[[group]]==base_group,][, 
                                                c("bug", "med", "iqr", "n")]
    meds_plot$bug_ord <- with(meds_plot, 
                              reorder(bug, med, function(x) { x }, order=T))
    meds_plot$type <- paste("Connecting medians of\ntaxa from", base_group)

    ggplot2::update_geom_defaults("point", list(colour=NULL))
    p1 <- ggplot2::ggplot(data=meds_plot, ggplot2::aes(x=bug_ord, y=med)) +
            ggplot2::geom_point()
    p1 <- p1 + ggplot2::geom_boxplot(data=d, ggplot2::aes(x=bug, y=value,
                                                          color=changed),
                                     notch=notch) + color_scale + 
        ggplot2::facet_grid(as.formula(paste("~", group)), labeller=labeller)
    # To change dots and line colors here:
    p1 <- p1 + ggplot2::geom_point(ggplot2::aes(shape=type), color="#999999") +
            ggplot2::geom_line(ggplot2::aes(group=1, linetype=type),
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
