#include <Rcpp.h>
#include <pbbam/DataSet.h>

using namespace Rcpp;
using namespace PacBio::BAM;



std::string getFileNameFromDataset(std::string dataset_name) {
  DataSet d(dataset_name);
  auto name = d.Name();
  auto id = d.ResourceId();
  auto tags = d.Tags();
  auto extResources = d.ExternalResources(); // <ns0:ExternalResources>
  auto extResource = extResources.Children().at(0); // <ns0:ExternalResource Name="Fasta lambdaNEB"
  auto fname = extResource.Attribute("ResourceId");
  return fname;
}

//' Get Reference Fasta File Path From Dataset
//'
//' This function parses a referenceset.xml file and assuming it follows a "Standard"
//' format returns the path to the file name inside the file.
//'
//' @param filename The BAM file name
//' @export
//' @examples getFastaFileNameFromDatasetFile("referenceset.xml")
// [[Rcpp::export]]
std::string getFastaFileNameFromDatasetFile(std::string dataset_name) {
  return getFileNameFromDataset(dataset_name);
}

//' Get Aligned BAM Name
//'
//' This function parses a alignmentset.xml file and assuming it follows a "Standard"
//' format returns the path to the file name inside the file.
//'
//' @param filename The BAM file name
//' @export
//' @examples getBAMNameFromDatasetFile("alignmentset.xml")
// [[Rcpp::export]]
std::string getBAMNameFromDatasetFile(std::string dataset_name) {
  return getFileNameFromDataset(dataset_name);
}

