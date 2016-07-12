#include <Rcpp.h>

#include <boost/algorithm/string.hpp>

#include <pbbam/BamRecord.h>
#include <pbbam/Config.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/Orientation.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/PbiIndexedBamReader.h>
#include <pbbam/ReadGroupInfo.h>

#include <pbbam/virtual/VirtualPolymeraseReader.h>
#include <pbbam/virtual/VirtualPolymeraseBamRecord.h>
#include <pbbam/virtual/VirtualRegion.h>
#include <pbbam/virtual/VirtualRegionType.h>

#include <fstream>

using namespace Rcpp;
using namespace PacBio::BAM;


bool FileExists(const std::string& path)
{
  std::ifstream ifile(path.c_str());
  return ifile.good();
}

char getRandomBase() {
  char bases[] = {'A', 'C', 'G', 'T'};
  int position = floor(4 * runif(1))[0];
  return bases[position];
}

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
    case 'N':
    case 'n':
        v[i] = 6;
        break;
    default:
        std::string msg = "Character was not an A, C, G, T, N or -. Check that the data you are loading "
                          "(BAM or Reference Fasta) does not contain weird characters.  Strange Character was: ";
        Rcpp::stop(msg + bp);
      }
  }
  v.attr("class") = "factor";
  v.attr("levels") = CharacterVector::create("A", "C", "G", "T", "-", "N");
  return v;
}


/* Convert two arrays of previous and current basepair into a
   dinucleotide context vector.
*/
IntegerVector createDinucleotideFactorFromSeqs(const IntegerVector& curBP,
                                              const IntegerVector& nextBP) {
  IntegerVector ctx(curBP.size());
  for(int i = 0; i < curBP.size(); i++ ) {
    auto nxt = nextBP[i];
    auto cur = curBP[i];
    auto offset = cur == nxt ? 0 : 4;
    // Expensive check that can probably be removed?
    if ( cur < 0 || cur > 4 || nxt < 0 || nxt > 4) {
        Rcpp::stop("Can't make context factor with basepairs outside of 0-4");
    }
    ctx[i] = offset + nxt;
  }
  ctx.attr("class") = "factor";
  ctx.attr("levels") = CharacterVector::create("AA", "CC", "GG", "TT", "NA", "NC", "NG", "NT");
  return ctx;
}


/* Given an aligned read and reference pair, we will generate
   breakpoints to divide the pairing into a smaller window everytime we
   see a window of at least newWindowSize
*/
std::vector<size_t> _findBreakPoints(const std::string& read,
                                     const std::string& ref,
                                     int newWindowSize) {

  if (read.length() != ref.length()) {
    Rcpp::stop("Read and reference had different lengths!");
  }
  std::vector<size_t> breaks;

  // Sample down if needed
  auto n = read.length();
  int lastBreak = -1;
  bool lastBaseMatch = false;
  // Start a new window everytime we find two bases that perfectly match every newWindowSize
  for(int i=0; i < n; i++) {
    bool curBaseMatch = (read[i] != '-' && ref[i] == read[i]) ? true : false;
    if (curBaseMatch && lastBaseMatch &&
        (lastBreak > 0 || // we are at the start
        (i - lastBreak) >= newWindowSize))
    {
      lastBreak = i;
      breaks.push_back(i);
      i += newWindowSize - 1;
    } else {
      lastBaseMatch = curBaseMatch;
    }
  }
  return breaks;
}

/* Given an aligned read and reference pair, we will subsample down to
   a smaller set of values
*/
std::pair<std::string, std::string> _sampleAndTrimSeqs(const std::string& read,
                                                       const std::string& ref,
                                                       int trimToLength) {

  if (read.length() != ref.length()) {
    Rcpp::stop("Read and reference had different lengths!");
  }
  // Sample down if needed
  auto n = read.length();
  size_t startRow, endRow;
  if (trimToLength !=0 && n > trimToLength) {
      auto maxRow = n - trimToLength + 1;
      startRow = static_cast<int>(maxRow * R::runif(0,1));
      endRow  = std::min(n, startRow + trimToLength);
  } else {
    startRow =0;
    endRow = n;
  }

  /* Let's make sure we don't start or end with a gap.
     We want the first and last two
     This could be made more efficient.
  */
  while (startRow < (n - 1) && (
         read[startRow] == '-' ||
         ref[startRow] == '-' ||
         ref[startRow + 1] == '-' ||
         read[startRow + 1] == '-') ) {
    startRow++;
  }
  while (endRow > 1 && (
         read[endRow] == '-' ||
         ref[endRow] == '-' ||
         read[endRow-1] == '-' ||
         ref[endRow-1] == '-')){
    endRow--;
  }

  if (read[endRow] == '-' || ref[endRow] == '-' ||  // we might have gone all the way to the end
      read[startRow] == '-' || ref[startRow] == '-' ||
      read[endRow-1] == '-' || ref[endRow-1] == '-' ||  // we might have gone all the way to the end
      read[startRow+1] == '-' || ref[startRow+1] == '-' ||
      endRow - startRow < 10) {
    Rcpp::stop("Trying to sample reads down wound up with \
                                  a sequence of size < 10 or an alignment that \
                                  started/ended with a gap at or immediately \
                                  adjacent to the end.  Are there a lot of \
                                  gaps or small sequences in this data?");
  }

  // Now let's subset
  int length = endRow - startRow + 1;
  auto new_read = read.substr(startRow, length);
  auto new_ref = ref.substr(startRow, length);
  return std::pair<std::string, std::string>(new_read, new_ref);
}



//' Load BAM header
//'
//' This function loads the header of a BAM file into a List.  It will attempt to load the
//' list of sequences, program used and read groups as well as the version of the
//' BAM file, skipping over any information that is not present.
//'
//' @param filename The BAM file name
//' @export
//' @examples loadHeader("~git/pbbam/tests/data/dataset/bam_mapping_1.bam")
// [[Rcpp::export]]
List loadHeader(std::string filename) {
  BamReader br(filename);
  auto head = br.Header();
  auto version = head.Version();
  auto df =  List::create( Named("version") = version);

  auto seqs = head.Sequences();
  if (seqs.size() >0) {
    auto n = seqs.size();
    CharacterVector names(n);
    CharacterVector lengths(n);
    for(int i=0; i < n; i++) {
      auto seq = seqs.at(i);
      names[i] = seq.Name();
      lengths[i] = seq.Length();
    }
    df["sequences"] = DataFrame::create(Named("reference") = names,
                                        Named("length") = lengths);
  }

  auto rgs = head.ReadGroups();
  if (rgs.size() > 0) {
    CharacterVector ids(rgs.size());
    CharacterVector movie(rgs.size());
    CharacterVector plat(rgs.size());
    CharacterVector chem(rgs.size());
    CharacterVector frameRate(rgs.size());
    CharacterVector binding(rgs.size());
    CharacterVector baseVersion(rgs.size());
    for(int i = 0; i < rgs.size(); i++) {
      auto rg = rgs.at(i);
      ids[i] = rg.Id();
      movie[i] = rg.MovieName();
      plat[i] = rg.Platform();
      frameRate[i] = rg.FrameRateHz();
      binding[i] = rg.BindingKit();
      baseVersion = rg.BasecallerVersion();
    }
    df["readgroups"] = DataFrame::create( _["movie"] = movie,
                                          _["platform"] = plat,
                                          _["chemistry"] = chem,
                                          _["framerate"] = frameRate,
                                          _["bindingkit"] = binding,
                                          _["basecaller"] = baseVersion);

  }


  auto comments = head.Comments();
  if (comments.size() > 0 ) {
    df["comments"] = wrap(comments);
  }

  auto programs = head.Programs();
  if (programs.size() > 0) {
    CharacterVector names(programs.size());
    CharacterVector version(programs.size());
    CharacterVector cmdLine(programs.size());
    for(int i=0; i < programs.size(); i++) {
      auto prog = programs.at(i);
      names[i] = prog.Name();
      version[i] = prog.Version();
      cmdLine[i] = prog.CommandLine();
    }
    df["programs"] = DataFrame::create(Named("name") = names,
                                       Named("version") = version,
                                       Named("cmdline") = cmdLine);
  }

  return df;
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
//' @param loadRQ Should we load the read quality? (Default = FALSE)
//' @export
//' @examples loadPBI("~git/pbbam/tests/data/dataset/bam_mapping_1.bam")
// [[Rcpp::export]]
DataFrame loadPBI(std::string filename,
                  bool loadSNR = false,
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
    IntegerVector mappedtId(wrap(mappedData.tId_));
    mappedtId = mappedtId + 1;
    try {
      BamReader br(filename);
      auto head = br.Header();
      if (head.Sequences().size() > 0) {
        mappedtId.attr("class") = "factor";
        mappedtId.attr("levels") = wrap(head.SequenceNames());
      }
    }catch(...) {}
    df["ref"] = mappedtId;
    df["tstart"] = mappedData.tStart_;
    df["tend"] = mappedData.tEnd_;
    df["astart"] = mappedData.aStart_;
    df["aend"] = mappedData.aEnd_;
    df["rc"] = LogicalVector(mappedData.revStrand_.begin(), mappedData.revStrand_.end());
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
//' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).
//' @param bamName The BAM file name to grab
//' @param indexedFastaName The name of the indexed fasta file this should come from.
//'
//' @return Returns a list of alignments as data frames.  If the IPD and Pulse Width are available, they will be columns in the returned data as well.
//' @export
// [[Rcpp::export]]
List loadDataAtOffsets(CharacterVector offsets, std::string bamName, std::string indexedFastaName) {
  try {

    if(!FileExists(indexedFastaName)) {
      stop("Fasta file does not exist or is not readable.");
    }

    if(!FileExists(bamName)) {
      stop("BAM file does not exist or is not readable.");
    }

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
          Rcpp::stop("Sequence and reference parts are of different size");
  // make this a list
        auto df = List::create(Named("read") = createFactorFromSeqString(seq),
                               Named("ref") = createFactorFromSeqString(ref));

        if (r.HasIPD()) {
            auto ipds = r.IPD(Orientation::NATIVE, true, true).Data();
            auto intarr = IntegerVector(ipds.size());
            for(int i = 0; i < ipds.size(); i++) {
              intarr[i] = ipds[i];
            }
            df["ipd"] = intarr;
        }

        if (r.HasPulseWidth()) {
            auto pws = r.PulseWidth(Orientation::NATIVE, true, true).Data();
            auto intarr = IntegerVector(pws.size());
            for(int i = 0; i < pws.size(); i++) {
              intarr[i] = pws[i];
            }
            df["pw"] = intarr;
        }

        if (r.HasSignalToNoise()) {
          auto snrs = r.SignalToNoise();
          df["snrA"] = rep(snrs[0], seq.size());
          df["snrC"] = rep(snrs[1], seq.size());
          df["snrG"] = rep(snrs[2], seq.size());
          df["snrT"] = rep(snrs[3], seq.size());
        }

        // Both are on hold pending bug fix #33721
        if(r.HasPkmid()) {
          df["pkmid"] = r.Pkmid(Orientation::NATIVE, true, true);
        }

        if (r.HasStartFrame()) {
          df["sf"] = r.StartFrame(Orientation::NATIVE, true, true);
        }

        df.attr("class") = "data.frame";
        df.attr("row.names") = seq_len(seq.size());
        results[i] = df;
        continue;
      } else{
        Rcpp::stop("No BAM record found at the given offset");
      }
    }
    return results;
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  }
  return List::create(Named("test") = 2);
}



//' Load BAM subreads as a list of data frames.
//'
//' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).
//' @param bamName The BAM file name to grab
//'
//' @return Returns a list of subreads as a list of character strings.
//' @export
// [[Rcpp::export]]
List loadSubreadsAtOffsets(CharacterVector offsets, std::string bamName) {

  try {
    if(!FileExists(bamName)) {
      stop("BAM file does not exist or is not readable.");
    }

    BamReader reader(bamName);
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
        std::string seq = r.Sequence();
        std::string name = r.FullName();
        auto fseq = List::create( Named("name") = name,
                                  Named("read") = seq);
        results[i] = fseq;
      } else{
        Rcpp::stop("No BAM record found at the given offset");
      }
    }
    return results;
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  }
  // Should never be called
  return List::create(Named("test") = 2);
}


//' Load a section of the reference genome as a character string.
//'
//' @param id The name of the sequence
//' @param start The start of the sequence
//' @param end The end of the sequence
//' @param indexedFastaName The name of the indexed fasta file this should come from.
//'
//' @return Returns a character vector of the sequence
//' @export
// [[Rcpp::export]]
CharacterVector loadReferenceWindow(std::string id, int start, int end, std::string indexedFastaName) {
  try {
    if(!FileExists(indexedFastaName)) {
      stop("Fasta file does not exist or is not readable.");
    }
    IndexedFastaReader fasta(indexedFastaName);
    // int -> int32_t (which is an int)
    Position s = start;
    Position e = end;
    return fasta.Subsequence(id, s, e);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  }
  return "failure";
}



//' Load BAM alignments as a list of list for the HMM model.
//'
//' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).
//' @param bamName The BAM file name to grab
//' @param indexedFastaName The name of the indexed fasta file this should come from.
//' @param trimToLength How much should we subsample the alignments?
//'
//' @return Returns a list of phase2datasets as data frames.
//' @export
// [[Rcpp::export]]
List loadHMMfromBAM(CharacterVector offsets,
                    std::string bamName,
                    std::string indexedFastaName,
                    int trimToLength = 140 ) {
  try {

    if (trimToLength < 0) {
      Rcout << "Trim to Length must be >= 0" << std::endl;
      return NULL;
    }

    IndexedFastaReader fasta(indexedFastaName);
    BamReader reader(bamName);

    // Always get reads in native orientation.
    auto orientation = Orientation::NATIVE;
    int n = offsets.size();
    List results(n);
    for(int i=0; i < n; i++) {
      // Check for interrupt
      if (i % 20 == 0 ) {
        Rcpp::checkUserInterrupt();
      }
      // back convert from string to long.
      std::string offset_string = as<std::string>(offsets[i]);
      long offset = std::stol(offset_string);
      reader.VirtualSeek(offset);
      BamRecord r;
      if (reader.GetNext(r)) {

        std::string read = r.Sequence(orientation, true, true);
        std::string ref = fasta.ReferenceSubsequence(r, orientation, true, true);
        // These should match and the first two and last two positions of the
        // alignment (no gaps at start or end).
        auto trimmed = _sampleAndTrimSeqs(read, ref, trimToLength);
        auto new_read = std::move(trimmed.first);
        auto new_ref = std::move(trimmed.second);

        // Trim out gaps
        boost::erase_all(new_ref, "-");
        boost::erase_all(new_read, "-");

        /* This is a bit brutal, but I think the nicer looking
           Rcpp sugar versions would involve more memory allocations.

           Note, since I am only using a dinucleotide context, I
           don't need to cut off the end, but I will just in
           case I want to make this trinucleotide later.
        */
        auto full_ref = createFactorFromSeqString(new_ref);
        size_t trimmed_size = full_ref.size() - 1;
        auto curBP = IntegerVector(trimmed_size);
        auto nextBP = IntegerVector(trimmed_size);
        for(int i = 0; i < trimmed_size; i++) {
          curBP[i] = full_ref[i];
          nextBP[i] = full_ref[i+1];
        }
        curBP.attr("class") = "factor";
        nextBP.attr("class") = "factor";
        curBP.attr("levels") = full_ref.attr("levels");
        nextBP.attr("levels") = full_ref.attr("levels");
        auto ctx = createDinucleotideFactorFromSeqs(curBP, nextBP);

        auto df =  List::create(Named("nextBP") = nextBP,
                                Named("currBP") = curBP,
                                Named("CTX") = ctx);

        // Now let's load the SNRs if applicable
        if (r.HasSignalToNoise()) {
            auto snrs = r.SignalToNoise();
            df["snrA"] = NumericVector(ctx.size(), snrs[0]);
            df["snrC"] = NumericVector(ctx.size(), snrs[1]);
            df["snrG"] = NumericVector(ctx.size(), snrs[2]);
            df["snrT"] = NumericVector(ctx.size(), snrs[3]);
          }
          df.attr("class") = "data.frame";
          df.attr("row.names") = seq_len(ctx.size());
          auto trimEnds = [](std::string val) {return val.substr(0, val.length() -1);};
          auto val = List::create(Named("covars") = df,
                                  Named("outcome") = CharacterVector(trimEnds(new_read)),
                                  Named("model") = CharacterVector(trimEnds(new_ref)));
        results[i] = val;
        continue;
      } else{
        Rcpp::stop("No BAM record found at the given offset");
      }
    }
    return results;
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  }
  return List::create(Named("test") = 2);
}


  List _generatePairFromString(const std::string& read, const std::string& ref, size_t brk1, size_t brk2) {
    // we should always have a break occur at a position with a matching base before, so we will grab that
    if (brk1 < 1) {
      Rcpp::stop("Break specified before 1! Shouldn't be possible as we require two matching bases at each breakpoint");
    }
    auto len = brk2 - brk1;
    auto new_read = read.substr(brk1, len);
    auto new_ref = ref.substr(brk1, len);
    // Trim out gaps
    boost::erase_all(new_ref, "-");
    boost::erase_all(new_read, "-");


    auto completeString = ref.substr(brk1, len + 1);
    boost::erase_all(completeString, "-");
    auto full_ref = createFactorFromSeqString(completeString);
    auto curBP = IntegerVector(full_ref.size() - 1);
    auto nextBP = IntegerVector(full_ref.size() - 1);
    for(int i = 0; i < (full_ref.size() - 1); i++) {
      curBP[i] = full_ref[i];
      nextBP[i] = full_ref[i + 1];
    }
    curBP.attr("class") = "factor";
    nextBP.attr("class") = "factor";
    curBP.attr("levels") = full_ref.attr("levels");
    nextBP.attr("levels") = full_ref.attr("levels");
    auto ctx = createDinucleotideFactorFromSeqs(curBP, nextBP);

    auto df =  List::create(Named("nextBP") = nextBP,
                            Named("currBP") = curBP,
                            Named("CTX") = ctx);

    df.attr("class") = "data.frame";
    df.attr("row.names") = seq_len(ctx.size());
    auto val = List::create(Named("covars") = df,
                            Named("outcome") = CharacterVector(new_read),
                            Named("model") = CharacterVector(new_ref));
    return val;
  }


//' Load BAM alignment from a single ZMW as a list of list for the HMM model.
//'
//' @param offsets The virtual file offsets to retrieve BAM records from (can be obtained from the index file based on loadPBI).
//' @param bamName The BAM file name to grab
//' @param indexedFastaName The name of the indexed fasta file this should come from.
//' @param windowBreakSize We generate a new "window" every time 2 basepairs are matching in a particular gap.
//' @param minSize What is the minimum window size necessary for return, if the window at the end is less than this, we drop it. (NOT CURRENTLY USED).
//' @return Returns a list of phase2datasets as data frames.
//' @export
// [[Rcpp::export]]
List loadSingleZmwHMMfromBAM(CharacterVector offsets,
                             std::string bamName,
                             std::string indexedFastaName,
                             int windowBreakSize = 140,
                             int minSize = 50 ) {
  try {

    if (windowBreakSize < minSize) {
      Rcout << "windowBreakSize must be >= minSize " << std::endl;
      return NULL;
    }

    IndexedFastaReader fasta(indexedFastaName);
    BamReader reader(bamName);

    // Always get reads in native orientation.
    auto orientation = Orientation::NATIVE;
    int n = offsets.size();
    if (n != 1) {
      Rcout << "You must specify exactly one offset" << std::endl;
    }
    // back convert from string to long.
    std::string offset_string = as<std::string>(offsets[0]);
    long offset = std::stol(offset_string);
    reader.VirtualSeek(offset);
    BamRecord r;
    if (reader.GetNext(r)) {
        std::string read = r.Sequence(orientation, true, true);
        std::string ref = fasta.ReferenceSubsequence(r, orientation, true, true);
        auto breaks = _findBreakPoints(read, ref, windowBreakSize);
        if (breaks.size() ==0) {
          Rcpp::stop("Could not find any break points in read pair.");
        }
        else {
          List results(breaks.size() - 1);
          for(int i=0; i < breaks.size() -1; i++) {
            auto temp = _generatePairFromString(read, ref, breaks[i], breaks[i+1]);
            results[i] = temp;
          }
          return results;
        }
    }
    else
    {
        Rcpp::stop("No BAM record found at the given offset");
    }
    } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  }
  return NULL;

  }




std::string replaceSuffix(const std::string& s,
                          const std::string& oldSuffix,
                          const std::string& newSuffix)
{
    if (s.substr(s.length() - oldSuffix.length(), oldSuffix.length()) != oldSuffix)
        stop("Unexpected suffix!");
    return s.substr(0, s.length() - oldSuffix.length()) + newSuffix;
}


//' Load regions table
//'
//' Load a "regions table" for a subreads.bam (and optionally, a
//' corresponding `aligned_subreads.bam`).
//'
//' Presently, this needs to traverse the entire BAM files---the PBI
//' doesn't contain any "region type" information.  If this becomes
//' crucial we might ' want to consider including it in the PBI.
//'
//' @param subreadsBamName the .subreads.bam file name.  The corresponding .scraps.bam must be available in the same directory
//' @return a data frame with columns HoleNumber, RegionType, RegionStart, RegionEnd, indicating the base extent of regions of each type.
//' @export
// [[Rcpp::export]]
DataFrame loadRegionsTable(const std::string& subreadsBamName)
{
    CharacterVector regionTypeLevels = CharacterVector::create(
        "ADAPTER", "BARCODE", "FILTERED", "SUBREAD",
        "HQREGION", "LQREGION", "ALIGNMENT");

    const std::map<VirtualRegionType, int> regionTypeToLevel {
        { VirtualRegionType::ADAPTER  , 1 },
        { VirtualRegionType::BARCODE  , 2 },
        { VirtualRegionType::FILTERED , 3 },
        { VirtualRegionType::SUBREAD  , 4 },
        { VirtualRegionType::HQREGION , 5 },
        { VirtualRegionType::LQREGION , 6 } };

    const std::string scrapsBamName = replaceSuffix(subreadsBamName, ".subreads.bam", ".scraps.bam");

    VirtualPolymeraseReader vpr(subreadsBamName, scrapsBamName);

    std::vector<int> holeNumber;
    std::vector<int> regionType;
    std::vector<int> startBase;
    std::vector<int> endBase;

    while (vpr.HasNext())
    {
        VirtualPolymeraseBamRecord record = vpr.Next();
        const auto virtualRegionsMap = record.VirtualRegionsMap();
        for (auto it : virtualRegionsMap)
        {
            const VirtualRegionType& rtype = it.first;
            const std::vector<VirtualRegion>& regions = it.second;
            for (const auto& region: regions) {
                holeNumber.push_back(record.HoleNumber());
                regionType.push_back(regionTypeToLevel.at(rtype));
                startBase.push_back(region.beginPos);
                endBase.push_back(region.endPos);
            }
        }
    }

    IntegerVector regionTypeR(regionType.begin(), regionType.end());
    regionTypeR.attr("class") = "factor";
    regionTypeR.attr("levels") = regionTypeLevels;

    DataFrame rv = DataFrame::create(
        Named("HoleNumber")  = holeNumber,
        Named("RegionType")  = regionTypeR,
        Named("RegionStart") = startBase,
        Named("RegionEnd")   = endBase);
    return rv;
}
