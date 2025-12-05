// This program takes the outfile from the flight path out toward the
// neutrino gravitational focus and creates a stacked histogram 



void histoCombiner () {  // file for opening
  
  TFile* myFile = TFile::Open("100000Orbits_100.000000kg_300.000000secondTimeStep_EllipticalOrbit_Closest=15.000000_Furthest=850000.root");
  //T33.0169File* myFile = TFile::Open("10000Orbits_100.000000kg_1800.000000secondTimeStep_FileOrbit.root");

  double flightTime = 38.8957;
  double scaleToYear = 365.24/flightTime;

  double scalar = 1e-5*scaleToYear;
  
  const double AUperRSolar = 1/215.032;
  const double maxR = 35;

  // Holder Histos
  TH1D* nRadHist = new TH1D("nRadHist","",218,-0.5,216.5);
  TH1D* solarRadHist = new TH1D("solarRadHist","",218,-0.5,216.5);
  TH1D* cosmicRadHist = new TH1D("cosmicRadHist","",218,-0.5,216.5);
  TH1D* RadHist = new TH1D("RadHist","",218,-0.5,216.5);
  
  TH1D* nArbRadHist = new TH1D("nArbRadHist","",maxR + 2,-0.5,maxR+0.5);
  TH1D* solarArbRadHist = new TH1D("solarArbRadHist","",maxR + 2,-0.5,maxR+0.5);
  TH1D* cosmicArbRadHist = new TH1D("cosmicArbRadHist","",maxR + 2,-0.5,maxR+0.5);

  TNtuple* neutrinoTuple = new TNtuple();

  /*
  myFile->GetObject("neutrinoRadiusHistogram",nRadHist);
  myFile->GetObject("neutrinoOnlyRadiusHistogram",nOnlyRadHist);
  myFile->GetObject("sub35NeutrinoRadiusHistogram",n35RadHist);
  myFile->GetObject("sub35NeutrinoOnlyRadiusHistogram",n35OnlyRadHist);
  myFile->GetObject("neutrinoFocusHistogram",nFocusRadHist);
  myFile->GetObject("neutrinoOnlyFocusHistogram",nOnlyFocusRadHist);
  myFile->GetObject("solarBackRadiusHistogram",solarRadHist);
  myFile->GetObject("sub35SolarBackRadiusHistogram",solar35RadHist);
  myFile->GetObject("cosmicBackRadiusHistogram",cosmicRadHist);
  myFile->GetObject("sub35CosmicBackRadiusHistogram",cosmic35RadHist);
  myFile->GetObject("radiusHistogram",RadHist);
  */
  
  myFile->GetObject("myTuple",neutrinoTuple);
  
  
  //int solarEntries = solarRadHist -> GetEntries(); 
  //int cosmicEntries = cosmicRadHist -> GetEntries(); 


  TCanvas *c0 = new TCanvas("c0","NeutrinoHistogram",1920,1080);
  c0->Divide(2,2);
  c0->cd(1);
  //gStyle -> SetOptStat(001000000);
  gStyle -> SetOptStat("i");
  //c0 -> SetLogy();

  neutrinoTuple -> Draw("neutrinoRadius","neutrinoRadius > 1 && neutrinoRadius < 35");
  auto htemp0 = (TH1F*)gPad->GetPrimitive("htemp");
  htemp0 -> Scale(scalar);
  //htemp0 -> SetTitle("Neutrino Count for Venus to 3 RSol Orbit;Distance From Sun (Solar Radii);Neutrinos/100kg/bin/orbit");
  //htemp0 -> SetTitle("Neutrino Count for ACO Orbit;Distance From Sun (Solar Radii);Neutrinos/100kg/bin/orbit");
  htemp0 -> SetTitle("Neutrino Count for Mercury to 15 Solar Radii;Distance From Sun (Solar Radii);Neutrinos/100kg/bin/year");
  htemp0 -> Draw("HIST");
  gPad->Update();
  TPaveStats *st0 = (TPaveStats*)htemp0->FindObject("stats");
  st0->SetX1NDC(0.70);
  st0->SetX2NDC(0.85);
  st0->SetY1NDC(0.80);
  st0->SetY2NDC(0.85);
  //c0 -> Print("nArbRadHist.svg");

  //TCanvas *c1 = new TCanvas("c1","CosmicHistogram",900,1200);
  //c1 -> SetLogy();

  c0->cd(2);
  neutrinoTuple -> Draw("cosmicBackRadius","cosmicBackRadius > 1 && cosmicBackRadius < 35");
  auto htemp1 = (TH1F*)gPad->GetPrimitive("htemp");
  htemp1 -> Scale(scalar);
  htemp1 -> SetTitle("Cosmic Ray Background for Mercury to 15 Solar Radii;Distance From Sun (Solar Radii);Cosmic Rays/bin/year");
  //htemp1 -> SetTitle("Cosmic Ray Background for ACO Orbit;Distance From Sun (Solar Radii);Neutrinos/100kg/bin/orbit");
  htemp1 -> Draw("HIST");
  gPad->Update();
  TPaveStats *st1 = (TPaveStats*)htemp1->FindObject("stats");
  st1->SetX1NDC(0.70);
  st1->SetX2NDC(0.85);
  st1->SetY1NDC(0.80);
  st1->SetY2NDC(0.85);
  //c1 -> Print("nArbRadHist.svg");

  //TCanvas *c2 = new TCanvas("c2","SolarHistogram",900,1200);
  //c2 -> SetLogy();
  c0->cd(3);

  neutrinoTuple -> Draw("solarBackRadius","solarBackRadius > 1 && solarBackRadius < 35");
  auto htemp2 = (TH1F*)gPad->GetPrimitive("htemp");
  htemp2 -> Scale(scalar);
  htemp2 -> SetTitle("Solar Wind Background for Mercury to 15 Solar Radii;Distance From Sun (Solar Radii);Solar Wind Protons/bin/year");
  //htemp2 -> SetTitle("Solar Wind Background for ACO Orbit;Distance From Sun (Solar Radii);Neutrinos/100kg/bin/orbit");
  htemp2 -> Draw("HIST");
  gPad->Update();
  TPaveStats *st2 = (TPaveStats*)htemp2->FindObject("stats");
  st2->SetX1NDC(0.70);
  st2->SetX2NDC(0.85);
  st2->SetY1NDC(0.80);
  st2->SetY2NDC(0.85);
  //c2 -> Print("nArbRadHist.png");
  

  cout << "I put the histograms in!.\n";


  
  





  
  THStack* myStack = new THStack("myStack","Stacked Events per Year;Radius(Solar Radii);Events/Solar Radii/year");

  TLegend *EnerLegend = new TLegend(0.45,0.65,0.85,0.85);
  EnerLegend -> SetHeader("Legend");

  
  htemp0->SetFillColor(kRed);
  htemp0 ->Rebin(2);
  myStack->Add(htemp0);
  
  htemp1->SetFillColor(kYellow);
  htemp1->Rebin(2);
  myStack->Add(htemp1);
  
  htemp2->SetFillColor(kBlue);
  htemp2->Rebin(2);
  myStack -> Add(htemp2);

  
  //cosmicRadHist->SetFillColor(kRed);
  //cosmicRadHist->Rebin(2);
  //myStack->Add(cosmicRadHist);
  
  //solarRadHist ->SetFillColor(kYellow);
  //solarRadHist->Rebin(2);
  //myStack->Add(solarRadHist);
  
  //nRadHist->SetFillColor(kBlue);
  //nRadHist->Rebin(2);
  //myStack -> Add( nRadHist );


  //TCanvas *c3 = new TCanvas("c3","Stacked",700,900);
  c0->cd(4);
  //c3 -> SetLogy();
  myStack->Draw("HIST");
  //myStack -> Scale(1/10000);
  //  myStack->GetXaxis()->
  
  EnerLegend -> AddEntry(nRadHist,"Neutrino Signal");
  EnerLegend -> AddEntry(solarRadHist,"Solar Background");
  EnerLegend -> AddEntry(cosmicRadHist,"Cosmic Background");
  EnerLegend -> Draw();




  /*
  THStack* myStack1 = new THStack("myStack1","Events;Radius(Solar Radii);Events/Bin");

  TLegend *EnerLegend1 = new TLegend(0.45,0.65,0.85,0.85);
  EnerLegend1 -> SetHeader("Legend");
  
  cosmic35RadHist->SetFillColor(kRed);
  //cosmicRadHist->Rebin(2);
  myStack1->Add(cosmic35RadHist);
  
  solar35RadHist ->SetFillColor(kYellow);
  //solarRadHist->Rebin(2);
  myStack1->Add(solar35RadHist);
  
  n35RadHist->SetFillColor(kBlue);
  //n35RadHist->Rebin(2);
  myStack1-> Add( n35RadHist );




  TCanvas *c2 = new TCanvas("c2","Stacked",700,900);
  myStack1->Draw("");

  EnerLegend1 -> AddEntry(n35RadHist,"Neutrino Signal");
  EnerLegend1 -> AddEntry(solar35RadHist,"Solar Background");
  EnerLegend1 -> AddEntry(cosmic35RadHist,"Cosmic Background");
  EnerLegend1 -> Draw();
  */



  c0 -> Print("Histograms.pdf");






}











