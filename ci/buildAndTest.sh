#!/bin/bash
# This script will be run within the Docker container to build and test.
R CMD INSTALL . && Rscript -e 'library(testthat); test_file("tests/testthat.R");'
