#include <Rcpp.h>

#include <pbbam/BamRecord.h>
#include <pbbam/Config.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/Orientation.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/PbiIndexedBamReader.h>
#include <pbbam/ReadGroupInfo.h>

using namespace Rcpp;
using namespace PacBio::BAM;

/* Convert the string into a factor variable, which is really just an integer
   vector with assigned levels.
*/
IntegerVector createFactorFromSeqString(const std::string& seq) {
  IntegerVector v(seq.length());
  for(int i = 0; i < seq.length(); i++ ) {
    char bp = seq[i];
    switch (bp) {
      case 'A':
      case 'a':
        v[i] = 1;
        break;
      case 'C':
      case 'c':
        v[i] = 2;
        break;
      case 'G':
      case 'g':
        v[i] = 3;
        break;
      case 'T':
      case 't':
        v[i] = 4;
        break;
      case '-':
        v[i] = 5;
        break;
      default:
        throw new std::out_of_range("Character was not an A, C, G, T or -");
      }
  }
  v.attr("class") = "factor";
  v.attr("levels") = CharacterVector::create("A", "C", "G", "T", "-");
  return v;
}


//' Load PBI BAM index file
//'
//' This function loads a pbi index file into a dataframe.  Depending on the
//' number of attributes present, it will either load just the basic data, or optionally
//' the mapping and barcode data.
//'
//' The original BAM file can also be read to gather additional covariates such as the SNR, read quality 
//' and number of passes, though this may take longer.
//'
//' @param filename The BAM file name (without .pbi)
//' @param loadSNR Should we load the four channel SNR data? (Default = FALSE)
//' @param loadNumPasses Should we load the number of passes data? (Default = FALSE)
//' @param loadRQ Shouold we load the read quality? (Default = FALSE)
//' @export
//' @examples loadpbi("~git/pbbam/tests/data/dataset/bam_mapping_1.bam.pbi")
// [[Rcpp::export]]
DataFrame loadpbi(std::string filename,
                  bool loadSNR= false,
                  bool loadNumPasses = false,
                  bool loadRQ = false
                  ) {
  /* Because the pbi file has different values depending on what it contains,
   * I am planning to construct it as a list of vectors, rather than a
   * dataframe directly, and then update the attributes to convert to a data frame
   */
  //"/Users/nigel/git/pbbam/tests/data/dataset/bam_mapping_1.bam.pbi"
  PbiRawData raw(filename + ".pbi");
  auto basicData = raw.BasicData();

  /* Because R can't handle longs, the safest way we will deal with this
   * is to convert them to strings. The most efficient would be to cast them
   * to doubles which I may implement later, but that is so ugly I can't bring
   * myself to do it now.
   */
   auto long_offsets = basicData.fileOffset_;
   std::vector<std::string> offsets(raw.NumReads());
   for(int i=0; i< raw.NumReads(); i++) {
      offsets[i] = std::to_string(long_offsets[i]);
   }

  // R data frames are basically lists
  auto df =  List::create( Named("hole") = basicData.holeNumber_,
                            Named("qstart") = basicData.qStart_,
                            Named("qend") = basicData.qEnd_,
                            Named("qual") = basicData.readQual_,
                            Named("offset") = offsets,
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
    // Convert to int to avoid raw vector assignment.
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

  // Load metrics from within the BAM?
  if (loadSNR || loadNumPasses || loadRQ) {
    /* Derek said every file in the pbi should be in the index, and every item
    should be ordered in the BAM so I am going to parse directly. */
    NumericVector snrA, snrC, snrG, snrT, np, rq;
    if (loadSNR) {
      snrA = NumericVector(raw.NumReads());
      snrC = NumericVector(raw.NumReads());
      snrG = NumericVector(raw.NumReads());
      snrT = NumericVector(raw.NumReads());
      df["snrA"] = snrA; df["snrC"] = snrC;
      df["snrG"] = snrG; df["snrT"] = snrT;
    }
    if(loadNumPasses) {
      np = NumericVector(raw.NumReads());
      df["np"] = np;
    }
    if (loadRQ) {
      rq = NumericVector(raw.NumReads());
      df["rq"] = rq;
    }


    int i=0;
    BamReader br(filename);
    BamRecord r;
    while (br.GetNext(r) ) {      
      if (loadSNR)
      { 
          if (r.HasSignalToNoise()) {
            auto snrs = r.SignalToNoise();
            snrA[i] = snrs[0]; snrC[i] = snrs[1];
            snrG[i] = snrs[2]; snrT[i] = snrs[3];
          } else{
            auto na = NumericVector::get_na();
            snrA[i] = na; snrC[i] = na;
            snrG[i] = na; snrT[i] = na;
          }
      }

      if (loadNumPasses) {
        if(r.HasNumPasses()) {
          np[i] = r.NumPasses();
        } else {
          np[i] = NumericVector::get_na();
        }
      }

      if (loadRQ) {
        if (r.HasReadAccuracy()) {
          rq[i] = r.ReadAccuracy();
        } else {
          rq[i] = NumericVector::get_na();
        }
      }
      i++;
    }
  }

  // Now let's make the list a data.frame
  df.attr("class") = "data.frame";
  df.attr("row.names") = seq_len(raw.NumReads());
  df.attr("bam.file") = filename;
  return df;
}


//' Load BAM alignments as a list of data frames.
//'
//' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadpbi).
//' @param bamName The BAM file name to grab
//' @param indexedFastaName The name of the indexed fasta file this should come from.
//'
//' @return Returns a list of alignments as data frames.
//' @export
// [[Rcpp::export]]
List loadDataAtOffsets(CharacterVector offsets, std::string bamName, std::string indexedFastaName) {
  try {

    IndexedFastaReader fasta(indexedFastaName);
    BamReader reader(bamName);

    // Always get reads in native orientation.
    auto orientation = Orientation::NATIVE;
    int n = offsets.size();
    List results(n);
    for(int i=0; i < n; i++) {
      // Check for interrupt
      if (i % 50 == 0 ) {
        Rcpp::checkUserInterrupt();
      }
      // back convert from string to long.
      std::string offset_string = as<std::string>(offsets[i]);
      long offset = std::stol(offset_string);
      reader.VirtualSeek(offset);
      BamRecord r;
      if (reader.GetNext(r)) {

        std::string seq = r.Sequence(orientation, true, true);
        std::string ref = fasta.ReferenceSubsequence(r, orientation, true, true);

        if (seq.size() != ref.size())
          throw std::runtime_error("Sequence and reference parts are of different size");
#if SNR
        // Generate SNR vectors
        auto snrs = r.SignalToNoise();
        NumericVector A(snrs[0], seq.length());
        NumericVector C(snrs[1], seq.length());
        NumericVector G(snrs[2], seq.length());
        NumericVector T(snrs[3], seq.length());
#endif
        DataFrame df = DataFrame::create(Named("read") = createFactorFromSeqString(seq),
                                         Named("ref") = createFactorFromSeqString(ref));
                                        // Named("snrA") = A,
                                        // Named("snrC") = C,
                                        // Named("snrG") = G,
                                        // Named("snrT") = T);
        results[i] = df;
        continue;
      } else{
        throw new std::out_of_range("No BAM record found at the given offset");
      }
    }
    return results;
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  }
  return List::create(Named("test") = 2);
}



