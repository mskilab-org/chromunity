
context("Unit testing chromunity operations")

library(chromunity)
library(GenomicRanges)
library(testthat)
library(pbmcapply)
library(arrow)
library(igraph)
library(skitools)
library(plyr)
library(MASS)

example_gr.path = system.file("extdata", "example_gr.rds", package = 'chromunity')
window_gr.path = system.file("extdata", "window_gr.rds", package = 'chromunity')
tiles_gr.path = system.file("extdata", "tiles_gr.rds", package = 'chromunity')
chromunity_out.path = system.file("extdata", "chromunity_out.rds", package = 'chromunity')
annotate_out.path = system.file("extdata", "annotate_out.rds", package = 'chromunity')
example_dir.path = system.file("extdata", "example_gr.rds", package = 'chromunity')

## Tests

test_that("parquet2gr", {
    example_dir1 = system.file("extdata/", package = 'chromunity')
    this.gr = parquet2gr(example_dir1, mc.cores = 2)
    expect_identical(class(this.gr)[1], "GRanges")
    example_dir2 = system.file("extdata/example_gr.rds", package = 'chromunity')
    expect_error(parquet2gr(example_dir2))
})

test_that("chromunity", {
    example_gr = readRDS(example_gr.path)
    this.chrom = chromunity(concatemers = example_gr, resolution = 1e5, window.size = 5e6, mc.cores = 2)
    expect_identical(class(this.chrom)[1], "Chromunity")
})

test_that("annotate", {
    example_chrom = readRDS(chromunity_out.path)
    binsets = gr2dt(example_chrom$binsets)
    binsets[, c := .N, by = bid]
    binsets = dt2gr(binsets[c > 1])
    this.annotate = annotate(binsets = binsets, concatemers = example_chrom$concatemers, mc.cores = 2)
    expect_identical(class(this.annotate)[1], "data.table")
})

test_that("synergy", {
    example_annot = readRDS(annotate_out.path)
    example_chrom = readRDS(chromunity_out.path)
    binsets = gr2dt(example_chrom$binsets)
    binsets[, c := .N, by = bid]
    binsets = dt2gr(binsets[c > 1])
    this.synergy = synergy(binsets = binsets, concatemers = example_chrom$concatemers, annotated.binsets = example_annot, resolution = 1e5, maxit = 10, mc.cores = 2)
    expect_identical(class(this.synergy)[1], "data.table")
})


