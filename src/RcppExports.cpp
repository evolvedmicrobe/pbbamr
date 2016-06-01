// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getFastaFileNameFromDatasetFile
std::string getFastaFileNameFromDatasetFile(std::string dataset_name);
RcppExport SEXP pbbamr_getFastaFileNameFromDatasetFile(SEXP dataset_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type dataset_name(dataset_nameSEXP);
    __result = Rcpp::wrap(getFastaFileNameFromDatasetFile(dataset_name));
    return __result;
END_RCPP
}
// getBAMNameFromDatasetFile
std::string getBAMNameFromDatasetFile(std::string dataset_name);
RcppExport SEXP pbbamr_getBAMNameFromDatasetFile(SEXP dataset_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type dataset_name(dataset_nameSEXP);
    __result = Rcpp::wrap(getBAMNameFromDatasetFile(dataset_name));
    return __result;
END_RCPP
}
// loadHeader
List loadHeader(std::string filename);
RcppExport SEXP pbbamr_loadHeader(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    __result = Rcpp::wrap(loadHeader(filename));
    return __result;
END_RCPP
}
// loadPBI
DataFrame loadPBI(std::string filename, bool loadSNR, bool loadNumPasses, bool loadRQ);
RcppExport SEXP pbbamr_loadPBI(SEXP filenameSEXP, SEXP loadSNRSEXP, SEXP loadNumPassesSEXP, SEXP loadRQSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< bool >::type loadSNR(loadSNRSEXP);
    Rcpp::traits::input_parameter< bool >::type loadNumPasses(loadNumPassesSEXP);
    Rcpp::traits::input_parameter< bool >::type loadRQ(loadRQSEXP);
    __result = Rcpp::wrap(loadPBI(filename, loadSNR, loadNumPasses, loadRQ));
    return __result;
END_RCPP
}
// loadDataAtOffsets
List loadDataAtOffsets(CharacterVector offsets, std::string bamName, std::string indexedFastaName);
RcppExport SEXP pbbamr_loadDataAtOffsets(SEXP offsetsSEXP, SEXP bamNameSEXP, SEXP indexedFastaNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< CharacterVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< std::string >::type bamName(bamNameSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    __result = Rcpp::wrap(loadDataAtOffsets(offsets, bamName, indexedFastaName));
    return __result;
END_RCPP
}
// loadSubreadsAtOffsets
List loadSubreadsAtOffsets(CharacterVector offsets, std::string bamName);
RcppExport SEXP pbbamr_loadSubreadsAtOffsets(SEXP offsetsSEXP, SEXP bamNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< CharacterVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< std::string >::type bamName(bamNameSEXP);
    __result = Rcpp::wrap(loadSubreadsAtOffsets(offsets, bamName));
    return __result;
END_RCPP
}
// loadReferenceWindow
CharacterVector loadReferenceWindow(std::string id, int start, int end, std::string indexedFastaName);
RcppExport SEXP pbbamr_loadReferenceWindow(SEXP idSEXP, SEXP startSEXP, SEXP endSEXP, SEXP indexedFastaNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    __result = Rcpp::wrap(loadReferenceWindow(id, start, end, indexedFastaName));
    return __result;
END_RCPP
}
// loadHMMfromBAM
List loadHMMfromBAM(CharacterVector offsets, std::string bamName, std::string indexedFastaName, int trimToLength);
RcppExport SEXP pbbamr_loadHMMfromBAM(SEXP offsetsSEXP, SEXP bamNameSEXP, SEXP indexedFastaNameSEXP, SEXP trimToLengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< CharacterVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< std::string >::type bamName(bamNameSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    Rcpp::traits::input_parameter< int >::type trimToLength(trimToLengthSEXP);
    __result = Rcpp::wrap(loadHMMfromBAM(offsets, bamName, indexedFastaName, trimToLength));
    return __result;
END_RCPP
}
// loadSingleZmwHMMfromBAM
List loadSingleZmwHMMfromBAM(CharacterVector offsets, std::string bamName, std::string indexedFastaName, int windowBreakSize, int minSize);
RcppExport SEXP pbbamr_loadSingleZmwHMMfromBAM(SEXP offsetsSEXP, SEXP bamNameSEXP, SEXP indexedFastaNameSEXP, SEXP windowBreakSizeSEXP, SEXP minSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< CharacterVector >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< std::string >::type bamName(bamNameSEXP);
    Rcpp::traits::input_parameter< std::string >::type indexedFastaName(indexedFastaNameSEXP);
    Rcpp::traits::input_parameter< int >::type windowBreakSize(windowBreakSizeSEXP);
    Rcpp::traits::input_parameter< int >::type minSize(minSizeSEXP);
    __result = Rcpp::wrap(loadSingleZmwHMMfromBAM(offsets, bamName, indexedFastaName, windowBreakSize, minSize));
    return __result;
END_RCPP
}
// loadRegionsTable
DataFrame loadRegionsTable(const std::string& subreadsBamName);
RcppExport SEXP pbbamr_loadRegionsTable(SEXP subreadsBamNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string& >::type subreadsBamName(subreadsBamNameSEXP);
    __result = Rcpp::wrap(loadRegionsTable(subreadsBamName));
    return __result;
END_RCPP
}
