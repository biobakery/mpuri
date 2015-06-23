# class levels for process_bug_names
class_levels <- c("s", "g", "f", "o", "c", "p", "k")

# not sure how to properly export this data 
# change colors here for coloring no change/increase/decrease
change_cols <- c("#0072B2", "#000000", "#D55E00")
change_cols <- c("#2C7BB6", "#FFFFBF", "#D7191C")
names(change_cols) <- c("decreased", "no change", "increased")
color_scale <- ggplot2::scale_colour_manual(name="change_cols",
                                            values=change_cols)

devtools::use_data(class_levels, color_scale, change_cols, internal = TRUE)
