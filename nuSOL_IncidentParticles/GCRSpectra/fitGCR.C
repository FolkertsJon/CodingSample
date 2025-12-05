// fitGCR.C — use TGraph::Eval, no helpers, relative-error weights,
// store E & Flux, integrate (midpoint), write total .mac + per-decade .macs,
// and ensure decade anchors 1,10,100,...,1e9 are explicitly evaluated.
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TKey.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TCollection.h" // TIter

static std::string outMacNameFromRoot(const char* rootname) {
  std::string s = rootname ? std::string(rootname) : std::string("output.root");
  size_t pos = s.rfind(".root");
  if (pos != std::string::npos) s.replace(pos, 5, ".mac");
  else s += ".mac";
  return s;
}

static std::string rangedName(const std::string& base, double lo, double hi) {
  size_t pos = base.rfind(".mac");
  char buf[128];
  std::snprintf(buf, sizeof(buf), "_%.0f-%.0f", lo, hi);
  if (pos != std::string::npos) return base.substr(0,pos) + std::string(buf) + base.substr(pos);
  return base + std::string(buf);
}

void fitGCR(const char* filename = "GCRProtons.root",
            double evalDiv = 1000.0) // 1000 turns GeV into MeV
{
  TFile* f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: could not open file '" << filename << "'\n";
    return;
  }

  std::vector<TGraphAsymmErrors*> graphs;
  for (TIter it(f->GetListOfKeys()); TKey* key = (TKey*)it(); ) {
    if (TString(key->GetClassName()) != "TGraphAsymmErrors") continue;
    TObject* obj = key->ReadObj();
    if (!obj) continue;
    if (!TString(obj->GetName()).EndsWith("_errtot")) continue;
    TGraphAsymmErrors* g = dynamic_cast<TGraphAsymmErrors*>(obj);
    if (g && g->GetN() > 0) graphs.push_back(g);
  }
  f->Close();
  if (graphs.empty()) {
    std::cerr << "Error: no *_errtot TGraphAsymmErrors found.\n";
    return;
  }

  std::vector<double> E;     // energies (MeV)
  std::vector<double> Flux;  // weighted averages scaled by evalDiv

  // We’ll step by ×1.1, but also inject exact decade anchors between steps.
  
  const double Emax = 1.0e9;
  double nextAnchor = 1.0e-2; // 1, 10, 100, ... up to 1e9

  for (double En = 1.0e-3; En < Emax; En *= 1.1) { // These are in GeV
    // 1) Evaluate at the current En
    {
      double num = 0.0, den = 0.0;
      for (size_t ig = 0; ig < graphs.size(); ++ig) {
        TGraphAsymmErrors* g = graphs[ig];
        g->Sort();

        const int n = g->GetN();
        if (n <= 0) continue;
	
        double xmin = +std::numeric_limits<double>::infinity();
        double xmax = -std::numeric_limits<double>::infinity();
        int idxNearest = -1;
        double bestAbs = +std::numeric_limits<double>::infinity();

        for (int i = 0; i < n; ++i) {
          double xi, yi;
          g->GetPoint(i, xi, yi);
          if (!std::isfinite(xi) || !std::isfinite(yi)) continue;
          if (xi < xmin) xmin = xi;
          if (xi > xmax) xmax = xi;
          double d = std::fabs(xi - En);
          if (d < bestAbs) { bestAbs = d; idxNearest = i; }
        }

        const double xEval = En;
        if (!(xmin < xmax)) continue; // skip if graph is funny
        if (xEval < xmin || xEval > xmax) continue; // skip if energy is outside range of graph
        if (idxNearest < 0) continue; // skip if no nearest index found

        double xi_near, yi_near;
        g->GetPoint(idxNearest, xi_near, yi_near);
        if (!(xi_near > 0.0) || !std::isfinite(yi_near) ) continue;


	// calculate weight for weighted average using average of square of relative errors
        const double exl = g->GetErrorXlow(idxNearest);
        const double exh = g->GetErrorXhigh(idxNearest);
        const double eyl = g->GetErrorYlow(idxNearest);
        const double eyh = g->GetErrorYhigh(idxNearest);
        const double rel_exl = exl / xi_near;
        const double rel_exh = exh / xi_near;
        const double rel_eyl = eyl / yi_near;
        const double rel_eyh = eyh / yi_near;
        double w = (rel_exl*rel_exl + rel_exh*rel_exh + rel_eyl*rel_eyl + rel_eyh*rel_eyh) / 4.0;
        if (!(w > 0.0) || !std::isfinite(w)) continue;
	w = 1.0 / w; // divide by weight
	
        double val = g->Eval(xEval); // linear interpolation
        if (!std::isfinite(val)) continue;

        num += val * w;
        den += w;
      }
      if (den > 0.0) {
        const double avgScaled = (num / den) / evalDiv;
        if (E.empty() || En != E.back()) { // avoid duplicates
          E.push_back(En*evalDiv); // convert to MeV
          Flux.push_back(avgScaled); // already in per MeV
        }
      }
    }

    // 2) Inject decade anchors that lie in (En, En*1.1] (and also catch 1e9 at the end)
    double En_next = En * 1.1;
    // advance nextAnchor up to be >= En
    while (nextAnchor < En) nextAnchor *= 10.0;

    while (nextAnchor <= En_next && nextAnchor <= Emax) {
      double A = nextAnchor;
      // if we already stored En == A due to rounding, skip
      if (E.empty() || A != E.back()) {
        double numA = 0.0, denA = 0.0;
        for (size_t ig = 0; ig < graphs.size(); ++ig) {
          TGraphAsymmErrors* g = graphs[ig];
          g->Sort();

          const int n = g->GetN();
          if (n <= 0) continue;

          double xmin = +std::numeric_limits<double>::infinity();
          double xmax = -std::numeric_limits<double>::infinity();
          int idxNearest = -1;
          double bestAbs = +std::numeric_limits<double>::infinity();

          for (int i = 0; i < n; ++i) {
            double xi, yi;
            g->GetPoint(i, xi, yi);
            if (!std::isfinite(xi) || !std::isfinite(yi)) continue;
            if (xi < xmin) xmin = xi;
            if (xi > xmax) xmax = xi;
            double d = std::fabs(xi - A);
            if (d < bestAbs) { bestAbs = d; idxNearest = i; }
          }

          const double xEvalA = A;
          if (!(xmin < xmax)) continue;
          if (xEvalA < xmin || xEvalA > xmax) continue;
          if (idxNearest < 0) continue;

          double xi_near, yi_near;
          g->GetPoint(idxNearest, xi_near, yi_near);
          if (!(xi_near > 0.0) || !std::isfinite(yi_near) || yi_near == 0.0) continue;
	  
	// calculate weight for weighted average using average of square of relative errors
          const double exl = g->GetErrorXlow(idxNearest);
          const double exh = g->GetErrorXhigh(idxNearest);
          const double eyl = g->GetErrorYlow(idxNearest);
          const double eyh = g->GetErrorYhigh(idxNearest);
          const double rel_exl = exl / xi_near;
          const double rel_exh = exh / xi_near;
          const double rel_eyl = eyl / std::fabs(yi_near);
          const double rel_eyh = eyh / std::fabs(yi_near);
          double w = (rel_exl*rel_exl + rel_exh*rel_exh + rel_eyl*rel_eyl + rel_eyh*rel_eyh) / 4.0;
          if (!(w > 0.0) || !std::isfinite(w)) continue;
	  w = 1.0 / w; // divide by weight

          double valA = g->Eval(xEvalA);
          if (!std::isfinite(valA)) continue;

          numA += valA * w;
          denA += w;
        }
        if (denA > 0.0) {
          const double avgScaledA = (numA / denA) / evalDiv;
          E.push_back(A*evalDiv);
          Flux.push_back(avgScaledA);
        }
      }

      // move to the next decade anchor
      if (nextAnchor >= Emax) break;
      nextAnchor *= 10.0;
    }
  }

  // Ensure 1e9 is included if not already (in case the loop never stepped past it)
  if (E.empty() || E.back() < Emax) {
    double num = 0.0, den = 0.0;
    for (size_t ig = 0; ig < graphs.size(); ++ig) {
      TGraphAsymmErrors* g = graphs[ig];
      g->Sort();
      const int n = g->GetN();
      if (n <= 0) continue;

      double xmin = +std::numeric_limits<double>::infinity();
      double xmax = -std::numeric_limits<double>::infinity();
      int idxNearest = -1;
      double bestAbs = +std::numeric_limits<double>::infinity();

      for (int i = 0; i < n; ++i) {
        double xi, yi;
        g->GetPoint(i, xi, yi);
        if (!std::isfinite(xi) || !std::isfinite(yi)) continue;
        if (xi < xmin) xmin = xi;
        if (xi > xmax) xmax = xi;
        double d = std::fabs(xi - Emax);
        if (d < bestAbs) { bestAbs = d; idxNearest = i; }
      }

      const double xEval = Emax;
      if (!(xmin < xmax)) continue;
      if (xEval < xmin || xEval > xmax) continue;
      if (idxNearest < 0) continue;

      double xi_near, yi_near;
      g->GetPoint(idxNearest, xi_near, yi_near);
      if (!(xi_near > 0.0) || !std::isfinite(yi_near) || yi_near == 0.0) continue;

      const double exl = g->GetErrorXlow(idxNearest);
      const double exh = g->GetErrorXhigh(idxNearest);
      const double eyl = g->GetErrorYlow(idxNearest);
      const double eyh = g->GetErrorYhigh(idxNearest);
      const double rel_exl = exl / xi_near;
      const double rel_exh = exh / xi_near;
      const double rel_eyl = eyl / std::fabs(yi_near);
      const double rel_eyh = eyh / std::fabs(yi_near);
      double w = (rel_exl*rel_exl + rel_exh*rel_exh + rel_eyl*rel_eyl + rel_eyh*rel_eyh) / 4.0;
      if (!(w > 0.0) || !std::isfinite(w)) continue;
      w = 1.0 / w; // divide by weight

      double val = g->Eval(xEval);
      if (!std::isfinite(val)) continue;

      num += val * w;
      den +=  w;
    }
    if (den > 0.0) {
      const double avgScaled = (num / den) / evalDiv;
      E.push_back(Emax*evalDiv);
      Flux.push_back(avgScaled);
    }
  }

  if (E.size() < 2) {
    std::cerr << "Error: not enough sampled points to integrate.\n";
    return;
  }

  // Midpoint-rule integral over all E
  double totalFlux = 0.0;
  for (size_t i = 0; i + 1 < E.size(); ++i) {
    const double Ei = E[i], Ei1 = E[i+1];
    const double dE = Ei1 - Ei; if (!(dE > 0.0)) continue;
    const double avgFlux = (Flux[i] + Flux[i+1])/2.0;
    //const double Emid = 0.5 * (Ei + Ei1);
    //const double t = (Emid - Ei) / dE;
    //const double Fmid = Flux[i] + t * (Flux[i+1] - Flux[i]);
    totalFlux += avgFlux*dE;
  }

  // Write TOTAL file
  const std::string totalName = outMacNameFromRoot(filename);
  {
    std::ofstream ofs(totalName.c_str());
    if (!ofs) { std::cerr << "Error: could not open output file '" << totalName << "'\n"; return; }
    ofs.setf(std::ios::scientific);
    ofs.precision(8);

    ofs << "# Total flux is " << totalFlux << " m^-2 sr^-1 s^-1\n";
    ofs << "/gps/ene/type User\n";
    ofs << "/gps/hist/type energy\n";
    ofs << "/gps/hist/point 0 0\n";
    for (size_t i = 0; i < E.size(); ++i) ofs << "/gps/hist/point " << E[i] << " " << Flux[i] << "\n";
  }

  // Per-decade files
  std::vector<std::pair<double,double> > ranges;
  for (double lo = 1.0; lo < Emax; lo *= 10.0) ranges.push_back(std::make_pair(lo, lo*10.0));

  for (size_t r = 0; r < ranges.size(); ++r) {
    const double Rlo = ranges[r].first;
    const double Rhi = ranges[r].second;

    // integrate over overlap with [Rlo,Rhi]
    double subFlux = 0.0; bool anyOverlap = false;
    for (size_t i = 0; i + 1 < E.size(); ++i) {
      double a = E[i], b = E[i+1];
      if (b <= Rlo || a >= Rhi) continue;
      double seg_lo = (a > Rlo) ? a : Rlo;
      double seg_hi = (b < Rhi) ? b : Rhi;
      if (seg_hi > seg_lo) {
        anyOverlap = true;
        double dE = seg_hi - seg_lo;
        double Emid = 0.5 * (seg_lo + seg_hi);
        double t = (Emid - a) / (b - a);
        double Fmid = Flux[i] + t * (Flux[i+1] - Flux[i]);
        subFlux += Fmid * dE;
      }
    }
    if (!anyOverlap) continue;

    std::string rname = rangedName(totalName, Rlo, Rhi);
    std::ofstream ofs(rname.c_str());
    if (!ofs) { std::cerr << "Error: could not open output file '" << rname << "'\n"; continue; }
    ofs.setf(std::ios::scientific);
    ofs.precision(8);

    ofs << "# Total flux is " << subFlux << " m^-2 sr^-1 s^-1\n";
    ofs << "/gps/ene/type User\n";
    ofs << "/gps/hist/type energy\n";
    ofs << "/gps/hist/point 0 0\n";
    for (size_t i = 0; i < E.size(); ++i) {
      if (E[i] >= Rlo && E[i] <= Rhi) ofs << "/gps/hist/point " << E[i] << " " << Flux[i] << "\n";
    }
  }

  std::cout << "Wrote " << totalName << " and per-decade .mac files; "
            << E.size() << " points sampled (anchors included).\n";
}
