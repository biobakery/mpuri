context("Pretty print OTU names")

test_that("process_bug_name works for greengenes names", {
    bugname1 <- "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__Lactobacillus.s__"
    expect_equal(process_bug_name(bugname1, sep="\\."), 
                 "g__Lactobacillus.")

    bugname4 <- "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__.s__"
    expect_equal(process_bug_name(bugname4, sep="\\."), 
                 "f__Lactobacillaceae.")

    bugname5 <- "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__S24_C20.s__"
    expect_equal(process_bug_name(bugname5, sep="\\."), 
                 "f__Lactobacillaceae.g__S24_C20.")
})


test_that("process bug names works for greengenes names and other seps", {
    bugname1 <- "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__species"
    expect_equal(process_bug_name(bugname1, sep="; "), 
                 "g__Lactobacillus; s__species")

    bugname5 <- "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__S24_C20; s__"
    expect_equal(process_bug_name(bugname5, sep="; "), 
                 "f__Lactobacillaceae; g__S24_C20; ")
})


test_that("process_bug_name works with dashes", {
    bugname5 <- "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__S24-C20; s__"
    expect_equal(process_bug_name(bugname5, sep="; "), 
                 "f__Lactobacillaceae; g__S24-C20; ")
})


test_that("process_bug_name works with spaces", {
    bugname5 <- "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__S24 C20; s__"
    expect_equal(process_bug_name(bugname5, sep="; "), 
                 "f__Lactobacillaceae; g__S24 C20; ")
})


test_that("silva names work", {
    n <- "Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;uncultured organism_lib_9986"
    expect_equal(process_silva_name(n, sep=";"), 
                 "Streptococcus;uncultured organism_lib_9986;")
})
