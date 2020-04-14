
context("Unit testing chromunity operations")

library(chromunity)
library(GenomicRanges)

example_gr.path = system.file("extdata", "example_gr.rds", package = 'chromunity')
window_gr.path = system.file("extdata", "window_gr.rds", package = 'chromunity')
tiles_gr.path = system.file("extdata", "tiles_gr.rds", package = 'chromunity')
chromunity_out_gr.path = system.file("extdata", "chromunity_out_gr.rds", package = 'chromunity')

## Tests


test_that("chromunity", {
    example_gr = readRDS(example_gr.path)
    window_gr = readRDS(window_gr.path)
    tiles_gr = readRDS(tiles_gr.path)
    chromunity_out = chromunity(example_gr, which.gr = window_gr, tiles = tiles_gr, k.knn = 25, k.min = 3)
    expect_identical(class(chromunity_out)[1], "GRanges")
    expect_true(identical(colnames(values(chromunity_out)), c("read_idx", "tix", "count", "community", "num.memb")))
})


test_that("annotate_multimodal_communities", {
    chromunity_out_gr = readRDS(chromunity_out_gr.path)
    window_gr = readRDS(window_gr.path)
    annotate_c = annotate_multimodal_communities(chromunity_out_gr, which.gr = window_gr)
    expect_true(identical(colnames(values(annotate_c)), c("community", "read_idx", "tix", "count", "num.memb", "multimodal")))
})

