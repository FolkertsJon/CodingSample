void histMaker() {
  // Open the ROOT file
  TFile *file = TFile::Open("ROOTOutput/ACO.root", "READ");
    
  // Check if the file has been opened successfully
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file." << std::endl;
    return;
  }

  // Get the ntuple from the file
  TNtupleD *inputTuple = (TNtupleD*)file->Get("myTuple");
  TH1D *neutrinoHist = (TH1D*)file->Get("neutrinoWeighted");
    
  // Check if the ntuple has been retrieved successfully
  if (!inputTuple) {
    std::cerr << "Error retrieving myTuple." << std::endl;
    file->Close();
    return;
  }
  
  /*
    if (!neutrinoHist) {
    std::cerr << "Error retrieving neutrinoWeighted." << std::endl;
    file->Close();
    return;
  }
  */
  
  // Calculate the number of bins
  double binWidth = 1.0; // Increment size
  int nBins = static_cast<int>((217.5 - (-0.5)) / binWidth) + 1; // Adding 1 to include the upper edge

  // Create the histogram
  TH1D *hist0 = new TH1D("hist0", "Neutrino Measurement for ACO Orbit for Excited State ;Radius (R_{#odot}) ;Neutrino Measurement (mission^{-1} 100kg^{-1} R_{#odot}^{-1})", nBins, -0.5, 217.5);
  TH1D *hist1 = new TH1D("hist1", "Neutrino Measurement for ACO Orbit for Excited State ;Radius (R_{#odot}) ;Neutrino Rate (mission^{-1} 100kg^{-1} R_{#odot}^{-1})", 36, -0.5, 35.5);

  // Set the statistics options: only display the integral
  gStyle->SetOptStat("i");

  
  // Assuming your TNtupleD has a variable you want to histogramize called "varName"
  inputTuple->Draw("neutrinoRadius>>hist0", "neutrinoSignal");
  inputTuple->Draw("neutrinoRadius>>hist1", "neutrinoSignal");

  // Draw the histogram on a canvas
  TCanvas *c1 = new TCanvas("c1", "Canvas Title", 1920, 1080);
  c1->SetLogy();
  hist0->Draw();

  c1->SaveAs("ACOFar.pdf");
  hist1->Draw();

  c1->SaveAs("ACOClose.pdf");

  //std::cout << "Press Return to exit";
  //std::cin.ignore();
  
  // Clean up
  file->Close();
  //delete hist0;
  //delete hist1;
  //delete hist2;
  delete c1;
  delete file; // Make sure to delete the file object after closing it to prevent memory leaks.
}
