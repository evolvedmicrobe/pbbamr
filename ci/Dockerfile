FROM rocker/r-devel
USER docker
RUN Rscript -e 'install.packages("Rcpp")'
RUN Rscript -e 'install.packages("testthat")'
CMD "/bin/bash"
