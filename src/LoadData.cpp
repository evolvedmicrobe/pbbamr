#include <Rcpp.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/ReadGroupInfo.h>

using namespace Rcpp;
using namespace PacBio::BAM;

//' Load PBI BAM index file
//'
//' This function loads a pbi index file into a dataframe.  Depending on the
//' number of attributes present, it will either load just the basic data, or optionally
//' the mapping and barcode data.
//'
//' @param filename A character vector with the file name in it.
//'
//' @export
//' @examples loadpbi("~git/pbbam/tests/data/dataset/bam_mapping_1.bam.pbi")
// [[Rcpp::export]]
DataFrame loadpbi(std::string filename) {
  /* Because the pbi file has different values depending on what it contains,
   * I am planning to construct it as a list of vectors, rather than a
   * dataframe directly, and then update the attributes to convert to a data frame
   */
  //"/Users/nigel/git/pbbam/tests/data/dataset/bam_mapping_1.bam.pbi"
  PbiRawData raw(filename);
  auto basicData = raw.BasicData();
  // R data frames are basically lists
  auto df =  List::create( Named("hole") = basicData.holeNumber_,
                            Named("qstart") = basicData.qStart_,
                            Named("qend") = basicData.qEnd_,
                            Named("qual") = basicData.readQual_,
                            Named("offset") = basicData.fileOffset_,
                            // Convert to int to avoid having it become a "Raw" vector
                            Named("flag") =  std::vector<int>(basicData.ctxtFlag_.begin(), basicData.ctxtFlag_.end() ));
  // Add Mapping Data
  if (raw.HasMappedData()) {
    auto mappedData = raw.MappedData();

    df["ref"] = mappedData.tId_;
    df["tstart"] = mappedData.tStart_;
    df["tend"] = mappedData.tEnd_;
    df["astart"] = mappedData.aStart_;
    df["aend"] = mappedData.aEnd_;
    df["rc"] = mappedData.revStrand_;
    df["matches"] = mappedData.nM_;
    df["mismatches"] = mappedData.nMM_;
    // It can't handle int8_t
    df["mq"] = std::vector<int>(mappedData.mapQV_.begin(), mappedData.mapQV_.end());
    // Now get number of insertions and deletions
    auto insertions = IntegerVector(raw.NumReads());
    auto deletions = IntegerVector(raw.NumReads());
    for(int i=0; i < raw.NumReads(); i++) {
        auto cnts = mappedData.NumDeletedAndInsertedBasesAt(i);
        deletions[i] = cnts.first;
        insertions[i] = cnts.second;
    }
    df["inserts"] = insertions;
    df["dels"] = deletions;
  }
  // Add Barcode Data
  if(raw.HasBarcodeData()) {
    auto bcd = raw.BarcodeData();
    df["bcf"] = bcd.bcForward_;
    df["bcr"] = bcd.bcReverse_;
    df["bcqual"] = std::vector<int>(bcd.bcQual_.begin(), bcd.bcQual_.end());
  }
  // Now let's make the list a data.frame
  df.attr("class") = "data.frame";
  df.attr("row.names") = seq_len(raw.NumReads());
  df.attr("bam.file") = filename;
  return df;
}
