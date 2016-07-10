FROM rocker/r-devel
USER docker
RUN Rscript -e 'install.packages("Rcpp")'
RUN Rscript -e 'install.packages("testthat")'
RUN Rscript -e 'install.packages("data.table")'
CMD "/bin/bash"
