// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getFastaFileNameFromDatasetFile
std::string getFastaFileNameFromDatasetFile(std::string dataset_name);
RcppExport SEXP pbbamr_getFastaFileNameFromDatasetFile(SEXP dataset_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type dataset_name(dataset_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(getFastaFileNameFromDatasetFile(dataset_name));
    return rcpp_result_gen;
END_RCPP
}
// getBAMNamesFromDatasetFile
CharacterVector getBAMNamesFromDatasetFile(std::string dataset_name);
RcppExport SEXP pbbamr_getBAMNamesFromDatasetFile(SEXP dataset_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type dataset_name(dataset_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(getBAMNamesFromDatasetFile(dataset_name));
    return rcpp_result_gen;
END_RCPP
}
// loadHeader
List loadHeader(std::string filename);
RcppExport SEXP pbbamr_loadHeader(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(loadHeader(filename));
    return rcpp_result_gen;
END_RCPP
}
// loadPBI
DataFrame loadPBI(std::string filename, bool loadSNR, bool loadNumPasses, bool loadRQ, bool loadSC);
RcppExport SEXP pbbamr_loadPBI(SEXP filenameSEXP, SEXP loadSNRSEXP, SEXP loadNumPassesSEXP, SEXP loadRQSEXP, SEXP loadSCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< bool >::type loadSNR(loadSNRSEXP);
    Rcpp::traits::input_parameter< bool >::type loadNumPasses(loadNumPassesSEXP);
    Rcpp::traits::input_parameter< bool >::type loadRQ(loadRQSEXP);
    Rcpp::traits::input_parameter< bool >::type loadSC(loadSCSEXP);
    rcpp_result_gen = Rcpp::wrap(loadPBI(filename, loadSNR, loadNumPasses, loadRQ, loadSC));
    return rcpp_result_gen;
END_RCPP
}
// loadDataAtOffsets
List loadDataAtOffsets(CharacterVector offsets, std::string bamName, std::string indexedFastaName);
RcppExport SEXP pbbamr_loadDataAtOffsets(SEXP offsetsSEXP, SEXP bamNameSEXP, SEXP indexedFastaNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< std::string >::type bamName(bamNameSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    rcpp_result_gen = Rcpp::wrap(loadDataAtOffsets(offsets, bamName, indexedFastaName));
    return rcpp_result_gen;
END_RCPP
}
// loadSubreadsAtOffsets
List loadSubreadsAtOffsets(CharacterVector offsets, std::string bamName);
RcppExport SEXP pbbamr_loadSubreadsAtOffsets(SEXP offsetsSEXP, SEXP bamNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< std::string >::type bamName(bamNameSEXP);
    rcpp_result_gen = Rcpp::wrap(loadSubreadsAtOffsets(offsets, bamName));
    return rcpp_result_gen;
END_RCPP
}
// loadReferenceWindow
CharacterVector loadReferenceWindow(std::string id, int start, int end, std::string indexedFastaName);
RcppExport SEXP pbbamr_loadReferenceWindow(SEXP idSEXP, SEXP startSEXP, SEXP endSEXP, SEXP indexedFastaNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    rcpp_result_gen = Rcpp::wrap(loadReferenceWindow(id, start, end, indexedFastaName));
    return rcpp_result_gen;
END_RCPP
}
// loadHMMfromBAM
List loadHMMfromBAM(CharacterVector offsets, std::string bamName, std::string indexedFastaName, int trimToLength);
RcppExport SEXP pbbamr_loadHMMfromBAM(SEXP offsetsSEXP, SEXP bamNameSEXP, SEXP indexedFastaNameSEXP, SEXP trimToLengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< std::string >::type bamName(bamNameSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    Rcpp::traits::input_parameter< int >::type trimToLength(trimToLengthSEXP);
    rcpp_result_gen = Rcpp::wrap(loadHMMfromBAM(offsets, bamName, indexedFastaName, trimToLength));
    return rcpp_result_gen;
END_RCPP
}
// loadSingleZmwHMMfromBAM
List loadSingleZmwHMMfromBAM(CharacterVector offsets, std::string bamName, std::string indexedFastaName, int windowBreakSize, int minSize);
RcppExport SEXP pbbamr_loadSingleZmwHMMfromBAM(SEXP offsetsSEXP, SEXP bamNameSEXP, SEXP indexedFastaNameSEXP, SEXP windowBreakSizeSEXP, SEXP minSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< std::string >::type bamName(bamNameSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    Rcpp::traits::input_parameter< int >::type windowBreakSize(windowBreakSizeSEXP);
    Rcpp::traits::input_parameter< int >::type minSize(minSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(loadSingleZmwHMMfromBAM(offsets, bamName, indexedFastaName, windowBreakSize, minSize));
    return rcpp_result_gen;
END_RCPP
}
// loadRegionsTable
DataFrame loadRegionsTable(const std::string& subreadsBamName);
RcppExport SEXP pbbamr_loadRegionsTable(SEXP subreadsBamNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type subreadsBamName(subreadsBamNameSEXP);
    rcpp_result_gen = Rcpp::wrap(loadRegionsTable(subreadsBamName));
    return rcpp_result_gen;
END_RCPP
}
// loadExtras
DataFrame loadExtras(DataFrame& df, bool loadSNR, bool loadNumPasses, bool loadRQ, bool loadSC);
RcppExport SEXP pbbamr_loadExtras(SEXP dfSEXP, SEXP loadSNRSEXP, SEXP loadNumPassesSEXP, SEXP loadRQSEXP, SEXP loadSCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< bool >::type loadSNR(loadSNRSEXP);
    Rcpp::traits::input_parameter< bool >::type loadNumPasses(loadNumPassesSEXP);
    Rcpp::traits::input_parameter< bool >::type loadRQ(loadRQSEXP);
    Rcpp::traits::input_parameter< bool >::type loadSC(loadSCSEXP);
    rcpp_result_gen = Rcpp::wrap(loadExtras(df, loadSNR, loadNumPasses, loadRQ, loadSC));
    return rcpp_result_gen;
END_RCPP
}
// loadPBI2
DataFrame loadPBI2(std::string filename);
RcppExport SEXP pbbamr_loadPBI2(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(loadPBI2(filename));
    return rcpp_result_gen;
END_RCPP
}
// getReadReport
List getReadReport(std::string datasetname, std::string indexedFastaName);
RcppExport SEXP pbbamr_getReadReport(SEXP datasetnameSEXP, SEXP indexedFastaNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type datasetname(datasetnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    rcpp_result_gen = Rcpp::wrap(getReadReport(datasetname, indexedFastaName));
    return rcpp_result_gen;
END_RCPP
}
