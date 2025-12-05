
// Include necessary headers
#include "TSystem.h"
#include "TList.h"
#include "TFile.h"
#include "TParameter.h"
#include "TMath.h"
#include "TSystemDirectory.h"


void calcMeanAndStDev(){
  // Directory containing the .root files
  const char* dirname = "data";
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (!files) {
    cout << "No files found in the directory!" << endl;
    return;
  }

  // Variables for statistics
  int n = 0;
  double sum = 0, sumsq = 0;

  // Iterate over each file in the directory
  TSystemFile *file;
  TIter next(files);
  while ((file = (TSystemFile*)next())) {
    TString fname = file->GetName();
    if (!file->IsDirectory() && fname.EndsWith(".root")) {
      // Construct full file path
      TString fpath = TString::Format("%s/%s", dirname, fname.Data());

      // Open the ROOT file
      TFile* rootFile = TFile::Open(fpath);
      if (!rootFile || rootFile->IsZombie()) {
        cout << "Error opening file " << fpath << endl;
        continue;
      }

      // Get the TParameter object
      TParameter<double>* param = (TParameter<double>*)rootFile->Get("fractionalNuSeen");
      if (!param) {
        cout << "Parameter fractionalNuSeen not found in " << fpath << endl;
      } else {
        // Retrieve the value and add to the sums
        double val = param->GetVal();
        sum += val;
        sumsq += val * val;
        n++;
      }

      // Close the file
      rootFile->Close();
      delete rootFile;
    }
  }

  // Calculate the mean and sample standard deviation, if at least one file was processed
  if (n > 1) {
    double mean = sum / n;
    double variance = (sumsq - (sum * sum / n)) / (n - 1);
    double stddev = TMath::Sqrt(variance);
    double sem = stddev / TMath::Sqrt(n);

    // Output results
    cout << "Mean of fractionalNuSeen: " << mean << endl;
    cout << "Sample Standard Deviation of fractionalNuSeen: " << stddev << endl;    cout << "Standard Error of the Mean: " << sem << endl;
  } else {
    cout << "Insufficient data for statistical analysis." << endl;

  }
}
