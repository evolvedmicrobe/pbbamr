#include <fstream>
#include <memory>
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

#include "utils.h"

using namespace Rcpp;
using namespace PacBio::BAM;


class PbiFilteredIndex : public PbiRawData
{
private:

public:
    PbiFilteredIndex()
        : PbiRawData()
    {}

    PbiFilteredIndex(uint32_t numReads, PbiFile::Sections sections)
        : PbiRawData()
    {
        NumReads(numReads);
        FileSections(sections);
        // reserve space
        this->BasicData() = PbiRawBasicData(numReads);
        this->MappedData() = PbiRawMappedData(numReads);
        this->BarcodeData() = PbiRawBarcodeData(numReads);
    }

    // No conversions/coercions here---deferred until ToDataFrame
    void AddPbiRecord(const PbiRawData& idx, uint16_t bamFileNumber, size_t row)
    {
        auto& thisBd = this->BasicData();
        const auto& thatBd = idx.BasicData();
        thisBd.rgId_       .push_back(thatBd.rgId_[row]);
        thisBd.qStart_     .push_back(thatBd.qStart_[row]);
        thisBd.qEnd_       .push_back(thatBd.qEnd_[row]);
        thisBd.holeNumber_ .push_back(thatBd.holeNumber_[row]);
        thisBd.readQual_   .push_back(thatBd.readQual_[row]);
        thisBd.ctxtFlag_   .push_back(thatBd.ctxtFlag_[row]);
        thisBd.fileOffset_ .push_back(thatBd.fileOffset_[row]);
        thisBd.fileNumber_ .push_back(bamFileNumber);

        if (HasMappedData()) {
            auto& thisMd = this->MappedData();
            const auto& thatMd = idx.MappedData();
            thisMd.tId_       .push_back(thatMd.tId_[row]);
            thisMd.tStart_    .push_back(thatMd.tStart_[row]);
            thisMd.tEnd_      .push_back(thatMd.tEnd_[row]);
            thisMd.aStart_    .push_back(thatMd.aStart_[row]);
            thisMd.aEnd_      .push_back(thatMd.aEnd_[row]);
            thisMd.revStrand_ .push_back(thatMd.revStrand_[row]);
            thisMd.nM_        .push_back(thatMd.nM_[row]);
            thisMd.nMM_       .push_back(thatMd.nMM_[row]);
            thisMd.mapQV_     .push_back(thatMd.mapQV_[row]);
        }

        if (HasBarcodeData()) {
            auto& thisBcd = this->BarcodeData();
            const auto& thatBcd = idx.BarcodeData();
            thisBcd.bcForward_ .push_back(thatBcd.bcForward_[row]);
            thisBcd.bcReverse_ .push_back(thatBcd.bcReverse_[row]);
            thisBcd.bcQual_    .push_back(thatBcd.bcQual_[row]);
        }
    }

    // All
    DataFrame ToDataFrame(const std::vector<std::string>& bamFilenames,
                          const std::vector<std::string>& seqNames) const
    {
        // --- Basic info ---
        auto& basicData = this->BasicData();
        size_t numReads = basicData.fileOffset_.size();

        std::vector<uint16_t> fnraw = basicData.fileNumber_;
        auto fileNumber = IntegerVector(fnraw.cbegin(), fnraw.cend(),
                                        [](const uint16_t& val) { return static_cast<int>(val) + 1; });
        fileNumber.attr("class") = "factor";
        fileNumber.attr("levels") = bamFilenames;


        // Convert offsets to strings since we have no 64bit integers
        // in R.  Move to using doubles instead?  FIXME: this is
        // getting recognized as a FACTOR, we should just have it be a
        // CharacterVector
        auto& offsets = basicData.fileOffset_;
        std::vector<std::string> offsetsAsString(offsets.size());
        for (int i = 0; i < numReads; i++) {
            offsetsAsString[i] = std::to_string(offsets[i]);
        }

        // Convert flag to int to avoid having it become "Raw" vector
        auto flag = IntegerVector(basicData.ctxtFlag_.begin(),
                                  basicData.ctxtFlag_.end(),
                                  [](const uint8_t f){return static_cast<int>(f);});

        List df = List::create(
            Named("file")   = fileNumber,
            Named("hole")   = basicData.holeNumber_,
            Named("qstart") = basicData.qStart_,
            Named("qend")   = basicData.qEnd_,
            Named("qual")   = basicData.readQual_,
            Named("offset") = CharacterVector(offsetsAsString.begin(), offsetsAsString.end()),
            Named("flag")   = flag);


        // --- Mapped info ---
        if (HasMappedData()) {
            auto& mappedData = this->MappedData();
            IntegerVector mappedtId(wrap(mappedData.tId_));
            mappedtId = mappedtId + 1;
            df["ref"] = mappedtId;
            if (seqNames.size() > 0) {
                mappedtId.attr("class") = "factor";
                mappedtId.attr("levels") = wrap(seqNames);
            }

            df["tstart"] = IntegerVector(mappedData.tStart_.begin(), mappedData.tStart_.end());
            df["tend"] = IntegerVector(mappedData.tEnd_.begin(), mappedData.tEnd_.end());
            df["astart"] = IntegerVector(mappedData.aStart_.begin(), mappedData.aStart_.end());
            df["aend"] = IntegerVector(mappedData.aEnd_.begin(), mappedData.aEnd_.end());
            df["rc"] = LogicalVector(mappedData.revStrand_.begin(), mappedData.revStrand_.end());
            df["matches"] = mappedData.nM_;
            df["mismatches"] = mappedData.nMM_;
            // Convert to int to avoid raw vector assignment.
            df["mq"] = IntegerVector(mappedData.mapQV_.begin(), mappedData.mapQV_.end(),
                                     [](const uint8_t f){return static_cast<int>(f);});

            auto insertions = IntegerVector(numReads);
            auto deletions = IntegerVector(numReads);
            for(size_t i = 0; i < numReads; i++) {
                deletions[i] = mappedData.NumDeletedBasesAt(i);
                insertions[i] = mappedData.NumInsertedBasesAt(i);
            }
            df["inserts"] = insertions;
            df["dels"] = deletions;
        }

        if (HasBarcodeData()) {
            auto bcd = this->BarcodeData();
            df["bcf"] = bcd.bcForward_;
            df["bcr"] = bcd.bcReverse_;
            df["bcqual"] = std::vector<int>(bcd.bcQual_.begin(), bcd.bcQual_.end());
        }

        df.attr("class") = "data.frame";
        df.attr("row.names") = seq_len(numReads);

        return df;
    }
};


//' Load some extra data values that are not in the pbi file, only in
//' the BAM.  This method is intended to be used to fetch information
//' about a significantly reduced subset of the entire BAM file
//' contents.
//'
//' This function returns a dataframe with same nrow dimension as df, with columns containing the
//' requested data fields.
//'
//' @param df the result of a call to loadPBI, or a subset of the rows of such a result
//' @param loadSNR Should we load the four channel SNR data? (Default = FALSE)
//' @param loadNumPasses Should we load the number of passes data? (Default = FALSE)
//' @param loadRQ Should we load the read quality? (Default = FALSE)
//' @param loadSC Load the SC tag for a scraps.bam file? (Only possible if file ends with '.scraps.bam')
//' @export
// [[Rcpp::export]]
DataFrame loadExtras(DataFrame& df,
                     bool loadSNR = false,
                     bool loadNumPasses = false,
                     bool loadRQ = false,
                     bool loadSC = false)
{
    Environment base("package:base");
    Function levels = wrap(base["levels"]);
    CharacterVector fileNames = levels(df["file"]);
    size_t nrows = df.nrows();
    DataFrame extras;


    if (loadSC) {
        for (const auto& fn: fileNames) {
            if (!has_suffix(std::string(fn), ".scraps.bam"))
            {
                Rcpp::stop("Can only set loadSC = TRUE if the filename passed ends with .scraps.bam");
            }
        }
    }

    NumericVector snrA, snrC, snrG, snrT, np, rq;
    IntegerVector sc;
    if (loadSNR) {
        snrA = NumericVector(nrows);
        snrC = NumericVector(nrows);
        snrG = NumericVector(nrows);
        snrT = NumericVector(nrows);
        extras["snrA"] = snrA; extras["snrC"] = snrC;
        extras["snrG"] = snrG; extras["snrT"] = snrT;
    }
    if(loadNumPasses) {
        np = NumericVector(nrows);
        extras["np"] = np;
    }
    if (loadRQ) {
        rq = NumericVector(nrows);
        extras["rq"] = rq;
    }
    if (loadSC) {
        sc = IntegerVector(nrows);
        sc.attr("class") = "factor";
        sc.attr("levels") = CharacterVector::create("Adapter", "Barcode", "Filtered", "LQRegion", "HQRegion");
        extras["sc"] = sc;
    }

    // Prepare a BAM reader for each bam file referenced; we will then
    // seek within each BAM reader to get the records we need.
    std::vector<std::unique_ptr<BamReader>> readers;
    for (auto bf : fileNames) {
        readers.emplace_back(new BamReader(as<std::string>(bf)));
    }

    CharacterVector offsetsAsString = df["offset"];
    IntegerVector fileNumber = df["file"];
    IntegerVector holes = df["hole"];

    for (size_t i = 0; i < nrows; i++)
    {
        // Check for interrupt
        if (i % 50 == 0 ) {
            Rcpp::checkUserInterrupt();
        }

        std::string offsetString = as<std::string>(offsetsAsString[i]);
        long offset = std::stol(offsetString);
        int fileNo = (int)fileNumber[i] - 1;

        auto& reader = readers[fileNo];
        BamRecord r;

        reader->VirtualSeek(offset);
        if (!reader->GetNext(r)) {
            stop("Failure to seek to offset " + offsetString);
        }

        // Verify data matches, just check hole for now.
        // If these don't match, the index and order in the file may be different.
        if(r.HoleNumber() != holes[i]) {
            Rcpp::stop("Hole numbers did not match when parsing BAM to generate additional values for index.");
        }

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

        if (loadSC) {
            if(r.HasScrapRegionType()) {
                // levels are "Adapter", "Barcode", "Filtered", "LQRegion", "HQRegion");
                auto sc_type = r.ScrapRegionType();
                int sc_int = 0;
                switch (sc_type) {
                case VirtualRegionType::ADAPTER:
                    sc_int = 1;
                    break;
                case VirtualRegionType::BARCODE:
                    sc_int = 2;
                    break;
                case VirtualRegionType::FILTERED:
                    sc_int = 3;
                    break;
                case VirtualRegionType::LQREGION:
                    sc_int = 4;
                    break;
                case VirtualRegionType::HQREGION:
                    sc_int = 5;
                    break;
                default:
                    Rcpp::stop("Unknown sc tag found in BAM file while parsing.");
                }
                sc[i] = sc_int;
            } else{
                sc[i] = IntegerVector::get_na();
            }
        }
    }

    extras.attr("class") = "data.frame";
    extras.attr("row.names") = seq_len(nrows);
    return extras;
}


//' Load PBI BAM index file
//'
//' This function loads a pbi index file into a dataframe.  Depending on the
//' number of attributes present, it will either load just the basic data, or optionally
//' the mapping and barcode data.
//'
//' Input can be a BAM filename, or an XML dataset file---not a
//' FOFN.
//'
//' @export
// [[Rcpp::export]]
DataFrame loadPBI2(std::string filename)

{
    // Improvements over the original loadPBI:
    // - do filtering up front, and do it correctly---we couldn't
    //   handle rname filters previously
    // - relegate the gathering of data from the BAM file to
    //   loadExtras, which can be run on a small subset
    DataSet ds(filename);
    PbiFilteredIndex fi;

    if (has_suffix(filename, ".bam")) {
        // BAM.  No filtering.
        PbiRawData pbiIndex(ds);
        fi = PbiFilteredIndex(pbiIndex.NumReads(), pbiIndex.FileSections());
        for (size_t row = 0; row < pbiIndex.NumReads(); row++)
        {
            fi.AddPbiRecord(pbiIndex, 0, row);
        }

    } else if (has_suffix(filename, ".xml")) {
        // XML datasets are composites of multiple BAMs, and can
        // indicate filters, which we need to respect.

        // We expect this size annotation to be accurate, which will
        // only not be true if someone has manually goofed with the
        // XML file---not used API or dataset tool.
        uint32_t expectedReadCount = std::stoi(ds.Metadata().NumRecords());

        const PbiFilter pbiFilter = PbiFilter::FromDataSet(ds);

        const auto& bamFiles = ds.BamFiles();
        for (size_t i = 0; i < bamFiles.size(); i++) {
            const auto& bamFile = bamFiles[i];
            const auto pbiFilename = bamFile.PacBioIndexFilename();
            PbiRawData pbiIndex(pbiFilename);

            if (i == 0) { fi = PbiFilteredIndex(expectedReadCount, pbiIndex.FileSections()); }

            for (size_t row = 0; row < pbiIndex.NumReads(); row++) {
                if (pbiFilter.Accepts(pbiIndex, row)) {
                    fi.AddPbiRecord(pbiIndex, i, row);
                }
            }

        }

    } else {
        stop("Only .bam and .xml dataset input accepted by this call.");
    }

    // Need the seq names and the filenames in order to label things
    // as factors
    const auto& bamFiles = ds.BamFiles();
    std::vector<std::string> bamFilenames(bamFiles.size());
    std::transform(bamFiles.begin(), bamFiles.end(), bamFilenames.begin(),
                   [](const BamFile& bf) { return bf.Filename(); });

    BamReader brn(bamFilenames.front());
    auto rhead = brn.Header();
    auto seqNames = rhead.SequenceNames();

    // Check
    for (auto fn: bamFilenames) {
      BamReader br(fn);
      auto head = br.Header();
      if (head.Sequences().size() != seqNames.size()) {
        stop("BAM files do not have same number of references!");
      } else {
        auto sName = head.SequenceNames();
        if (!std::equal(seqNames.begin(), seqNames.end(), sName.begin())) {
          stop("BAM files have different reference names!");
        }
      }
    }

    DataFrame df =  fi.ToDataFrame(bamFilenames, seqNames);
    df.attr("bam.file") = filename;
    return df;
}
