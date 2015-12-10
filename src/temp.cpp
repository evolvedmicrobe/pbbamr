#include <Rcpp.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/ReadGroupInfo.h>


using namespace Rcpp;
using namespace PacBio::BAM;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
CharacterVector timesTwo(NumericVector x) {
  Rcpp::Rcout << "Armadillo version: "  << std::endl;
  DataSet ds2("/Users/nigel/pacbio/2xccs/m140830_011125_42129_c100699650030000001823144203261570_s1_p0.subreads.bam");
  EntireFileQuery query(ds2);
  PacBio::BAM::BamRecord rec;
  auto res = query.GetNext(rec);
  auto id = rec.FullName();
  Rcpp::Rcout << "Armadillo version: " << id << std::endl;
  return CharacterVector(res);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
