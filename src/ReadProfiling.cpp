// This file contains a series of methods designed to examine "per-read"
// metrics, those that require parsing of the BAM file as opposed to just
// examining the PBI files. In order to load each BAM into memory only once,
// we'll parse the file and pass it to each of the "error reporter" classes that
// derive from a base class that implements the general methods.

#include <memory>
#include <fstream>
#include <Rcpp.h>

#include <boost/algorithm/string.hpp>

#include <pbbam/BamRecord.h>
#include <pbbam/Config.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/Orientation.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/PbiIndexedBamReader.h>
#include <pbbam/ReadGroupInfo.h>


#include "utils.h"

using namespace Rcpp;
using namespace PacBio::BAM;

// This is the base class that all reporters derive from.
// to implement a new report just add a derived class and add it to the vector
// in `getReadReport`
class PerReadMetricReporter {
public:
  virtual void ConsumeRead(const BamRecord& r, const std::string ref, const std::string seq) = 0;
  virtual List ProduceReport() = 0;
  virtual std::string GetName() = 0;
};

// Collect counts of mismatches between bases.
class MismatchReport : public PerReadMetricReporter {
private:
  // Row is ref, column is read.
  int mismatches[BasePairs::N + 1][BasePairs::N + 1];

public:
  virtual void ConsumeRead(const BamRecord& r, const std::string ref, const std::string seq) {
    for(int i = 0; i < ref.size(); i++) {
      auto tbp = BPtoIndex(ref[i]);
      auto rbp = BPtoIndex(seq[i]);
      if (tbp != rbp && tbp != BasePairs::Gap  && rbp != BasePairs::Gap) {
          mismatches[tbp][rbp]++;
        }
    }
  }
  virtual List ProduceReport() {
      const auto N = BasePairs::N;
      auto arrSize = (N * N);
      IntegerVector refBases(arrSize);
      IntegerVector seqBases(arrSize);
      NumericVector cnts(arrSize);
      for(int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
          auto index = (i - 1)*N + (j - 1);
          refBases[index] = i;
          seqBases[index] = j;
          cnts[index] = mismatches[i][j];
        }
      }
      refBases.attr("class") = "factor";
      refBases.attr("levels") = CharacterVector::create("A", "C", "G", "T", "-", "N");
      seqBases.attr("class") = "factor";
      seqBases.attr("levels") = CharacterVector::create("A", "C", "G", "T", "-", "N");

      return DataFrame::create(_["ref"] = refBases,
                               _["read"] = seqBases,
                               _["cnts"] = cnts);
  }
  virtual std::string GetName() {
    return "Mismatches";
  }
};

class GapSizeReport : public PerReadMetricReporter {
private:
  // This is the largest gap size I will record, after this it's 50+
  const int maxGapSize = 50;
  std::vector<int> refGapCounts;
  std::vector<int> readGapCounts;
  void countGapSizes(const std::string& str, std::vector<int>& vec) {
    for (size_t i = 0; i < str.size(); i++) {
      if(str[i] == '-') {
        int gapSize = 0;
        while(i < str.size() && str[i] == '-') {
          gapSize++;
          i++;
        }
        // Only record the gap size if it isn't at the end, otherwise
        // we don't actually know what the size is.
        if (i != str.size() && (i - gapSize) > 0) {
          auto index = std::min(gapSize, maxGapSize);
          // 0 -> 1 encoding change
          vec[index - 1]++;
        }
      }
    }
  }
public:
  GapSizeReport() : refGapCounts(maxGapSize),
                    readGapCounts(maxGapSize) {};
  virtual void ConsumeRead(const BamRecord& r, const std::string ref, const std::string seq) {
    countGapSizes(ref, refGapCounts);
    countGapSizes(seq, readGapCounts);
  }
  virtual List ProduceReport() {
    return DataFrame::create(_["gapSize"] = seq_len(maxGapSize),
                             _["refCnts"] = refGapCounts,
                             _["readCnts"] = readGapCounts);
  }
  virtual std::string GetName() {
    return "gapSizes";
  }
};


//' Get Per Read Metrics
//'
//' This function loads a dataset and parses each read, passing it to a class
//' which collects metrics on each read.  It returns a list of data frames, one
//' for each metric analyzed.
//'
//' @param dataset The dataset/BAM file name.
//' @param indexedFastaName The fasta file used in the alignment.
//' @export
//' @examples
// [[Rcpp::export]]
List getReadReport(std::string datasetname, std::string indexedFastaName) {

  if(!FileExists(indexedFastaName)) {
    stop("Fasta file does not exist or is not readable.");
  }
  if(!FileExists(datasetname)) {
    stop("BAM file does not exist or is not readable.");
  }
  std::vector<std::unique_ptr<PerReadMetricReporter> > reporters;
  reporters.emplace_back(new MismatchReport());
  reporters.emplace_back(new GapSizeReport());
  IndexedFastaReader fasta(indexedFastaName);

  // Always get reads in native orientation.
  auto orientation = Orientation::NATIVE;

  DataSet ds(datasetname);
  const auto filter = PbiFilter::FromDataSet(ds);
  std::unique_ptr<PacBio::BAM::internal::IQuery> query(nullptr);
  if (filter.IsEmpty())
    query.reset(new EntireFileQuery(ds));
  else
    query.reset(new PbiFilterQuery(filter, ds));

  int interruptCounter = 0;
  for (const BamRecord& read : *query) {
    if (++interruptCounter % 100 == 0) Rcpp::checkUserInterrupt();
    std::string seq = read.Sequence(orientation, true, true);
    std::string ref = fasta.ReferenceSubsequence(read, orientation, true, true);
    if (seq.size() != ref.size()) Rcpp::stop("Sequence and reference are different sizes");
    for (auto& reporter : reporters) {
      reporter->ConsumeRead(read, ref, seq);
    }
  }

  auto toReturn = List();
  for(auto& reporter : reporters) {
    toReturn[reporter->GetName()] = reporter->ProduceReport();
  }
  return toReturn;
}
