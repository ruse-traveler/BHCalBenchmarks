// ----------------------------------------------------------------------------
// 'hcal_barrel_cluster_analysis_standalone.cxx'
// Derek Anderson
// 09.07.2023
//
// for generating ePIC BHCal benchmarks
// in a ROOT interactive session.
// ----------------------------------------------------------------------------

// c++ utilities
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>
// root classes
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
// dataframe related classes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>

// make common namespaces implicit
using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;

// set up aliases
using TH1Def = ROOT::RDF::TH1DModel;
using TH2Def = ROOT::RDF::TH2DModel;



// bhcal benchmarks -----------------------------------------------------------

int hcal_barrel_cluster_analysis_standalone(
  const string sInFile  = "./test/forTruthAssocTest.epic23090image.e10th45pip.podio.root",
  const string sOutFile = "test.root"
) {

  // turn on histogram errors
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  // open output file
  TFile* fOutput = new TFile(sOutFile.data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return -1;
  }

  // open dataframe
  RDataFrame frame("events", sInFile);

  // make sure file isn't empty
  auto events = frame.Count();
  if (events == 0) {
    cerr << "Error: No events found!" << endl;
    assert(events > 0);
  }

  // lambdas for analysis -----------------------------------------------------

  auto getParticleEnergy = [](const RVec<int> &types, const RVec<float> &energies) {
    float    energy = -1;
    uint64_t index  = 0;
    for (const int type : types) {
      if (type == 1) {
        energy = energies.at(index);
        break;
      }
      ++index;
    }
    return energy;
  };

  auto getLeadEnergy = [](const RVec<float> &energies) {
    float lead = -1;
    for (const float energy : energies) {
      if (energy > lead) {
        lead = energy;
      }
    }
    return lead;
  };

  auto getEnergySum = [](const RVec<float> &energies) {
    float sum = 0.;
    for (const float energy : energies) {
      sum += energy;
    }
    return sum;
  };

  auto getPercentDiffVec = [](const float a, const RVec<float> &b) {
    RVec<float> diff;
    for (uint64_t index = 0; index < b.size(); index++) {
      diff.push_back((a - b.at(index)) / a);
    }
    return diff;
  };

  auto getPercentDiff = [](const float a, const float b) {
    return (a - b) / b;
  };

  auto getSumFraction = [](const float a, const float b) {
    return a / (a + b);
  };

  // TODO template input vector
  auto getMultiplicity = [](const RVec<float> &collection) {
    return collection.size();
  };

  // define histograms --------------------------------------------------------

  // histogram titles
  const string eneTitle       = ";E [GeV];counts";
  const string diffTitle      = ";(E_{reco} - E_{par}) / E_{par};counts";
  const string multTitle      = ";multiplicity;counts";
  const string diffCvsPTitle  = ";E_{par} [GeV];E_{clust} [GeV];counts";
  const string diffLvsPTitle  = ";E_{par} [GeV];E_{clust}^{lead} [GeV];counts";
  const string diffSvsPTitle  = ";E_{par} [GeV];#SigmaE_{clust} [GeV];counts";
  const string sumVsFracTitle = ";#SigmaE_{BIC} / (#SigmaE_{BIC} + #SigmaE_{BHCal});#SigmaE_{BIC};counts";

  // histogram binning
  const tuple<int, double, double> eneBins  = make_tuple(50,  0.,   50.);
  const tuple<int, double, double> diffBins = make_tuple(100, -10., 10.);
  const tuple<int, double, double> fracBins = make_tuple(100, 0.,   10.);
  const tuple<int, double, double> multBins = make_tuple(100, 0.,   100.);

  // histogram definitions
  const vector<TH1Def> vecHistDefs1D = {
    TH1Def("hEnePar",         eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneHit",         eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneClust",       eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneLead",        eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneSum",         eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hDiffClust",      diffTitle.data(), get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hDiffLead",       diffTitle.data(), get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hDiffSum",        diffTitle.data(), get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hMultHit",        multTitle.data(), get<0>(multBins), get<1>(multBins), get<2>(multBins)),
    TH1Def("hMultClust",      multTitle.data(), get<0>(multBins), get<1>(multBins), get<2>(multBins)),
    TH1Def("hMultHitInClust", multTitle.data(), get<0>(multBins), get<2>(multBins), get<2>(multBins))
  };
  const vector<TH2Def> vecHistDefs2D = {
    TH2Def("hSumVsFrac", sumVsFracTitle.data(), get<0>(fracBins), get<1>(fracBins), get<2>(fracBins), get<0>(eneBins), get<1>(eneBins), get<2>(eneBins))
  };

  // run analysis -------------------------------------------------------------

  // define columns
  auto analysis = frame.Define("enePar",        getParticleEnergy, {"GeneratedParticles.type", "GeneratedParticles.energy"})
                       .Define("eneLeadHCal",   getLeadEnergy,     {"HcalBarrelClusters.energy"})
                       .Define("eneSumHCal",    getEnergySum,      {"HcalBarrelClusters.energy"})
                       .Define("eneSumECal",    getEnergySum,      {"EcalBarrelClusters.energy"})
                       .Define("diffClustHCal", getPercentDiffVec, {"enePar", "HcalBarrelClusters.energy"})
                       .Define("diffLeadHCal",  getPercentDiff,    {"enePar", "eneLeadHCal"})
                       .Define("diffSumHCal",   getPercentDiff,    {"enePar", "eneSumHCal"})
                       .Define("fracSumHxECal", getSumFraction,    {"eneSumECal", "eneSumHCal"})
                       .Define("multHit",       getMultiplicity,   {"HcalBarrelRecHits.energy"})
                       .Define("multClust",     getMultiplicity,   {"HcalBarrelClusters.energy"});

  // get 1D histograms
  //   TODO it might be nice to collect these into a vector
  //   and automate histogram operations...
  auto hEnePar         = analysis.Histo1D(vecHistDefs1D[0],  "enePar");
  auto hEneHit         = analysis.Histo1D(vecHistDefs1D[1],  "HcalBarrelRecHits.energy");
  auto hEneClust       = analysis.Histo1D(vecHistDefs1D[2],  "HcalBarrelClusters.energy");
  auto hEneLead        = analysis.Histo1D(vecHistDefs1D[3],  "eneLeadHCal");
  auto hEneSum         = analysis.Histo1D(vecHistDefs1D[4],  "eneSumHCal");
  auto hDiffClust      = analysis.Histo1D(vecHistDefs1D[5],  "diffClustHCal");
  auto hDiffLead       = analysis.Histo1D(vecHistDefs1D[6],  "diffLeadHCal");
  auto hDiffSum        = analysis.Histo1D(vecHistDefs1D[7],  "diffSumHCal");
  auto hMultHit        = analysis.Histo1D(vecHistDefs1D[8],  "multHit");
  auto hMultClust      = analysis.Histo1D(vecHistDefs1D[9],  "multClust");
  auto hMultHitInClust = analysis.Histo1D(vecHistDefs1D[10], "HcalBarrelClusters.nhits");

  // get 2D histograms
  auto hSumVsFrac = analysis.Histo2D(vecHistDefs2D[0], "fracSumHxECal", "eneSumECal");

  // normalize 1D histograms by integral
  const double intEnePar         = hEnePar         -> Integral();
  const double intEneHit         = hEneHit         -> Integral();
  const double intEneClust       = hEneClust       -> Integral();
  const double intEneLead        = hEneLead        -> Integral();
  const double intEneSum         = hEneSum         -> Integral();
  const double intDiffClust      = hDiffClust      -> Integral();
  const double intDiffLead       = hDiffLead       -> Integral();
  const double intDiffSum        = hDiffSum        -> Integral();
  const double intMultHit        = hMultHit        -> Integral();
  const double intMultClust      = hMultClust      -> Integral();
  const double intMultHitInClust = hMultHitInClust -> Integral();
  if(intEnePar         > 0.) hEnePar         -> Scale(1. / intEnePar);
  if(intEneHit         > 0.) hEneHit         -> Scale(1. / intEneHit);
  if(intEneClust       > 0.) hEneClust       -> Scale(1. / intEneClust);
  if(intEneLead        > 0.) hEneLead        -> Scale(1. / intEneLead);
  if(intEneSum         > 0.) hEneSum         -> Scale(1. / intEneSum);
  if(intDiffClust      > 0.) hDiffClust      -> Scale(1. / intDiffClust);
  if(intDiffLead       > 0.) hDiffLead       -> Scale(1. / intDiffLead);
  if(intDiffSum        > 0.) hDiffSum        -> Scale(1. / intDiffSum);
  if(intMultHit        > 0.) hMultHit        -> Scale(1. / intMultHit);
  if(intMultClust      > 0.) hMultClust      -> Scale(1. / intMultClust);
  if(intMultHitInClust > 0.) hMultHitInClust -> Scale(1. / intMultHitInClust);

  // make plots ---------------------------------------------------------------

  // define styles
  const tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> stylePar     = {1,   1, 2, 0, 24};
  const tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> styleHit     = {633, 1, 2, 0, 27};
  const tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> styleClust   = {417, 1, 2, 0, 25};
  const tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> styleLead    = {601, 1, 2, 0, 28};
  const tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> styleSum     = {617, 1, 2, 0, 46};
  const tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> styleFrac    = {1,   1, 2, 0, 1};
  const tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> styleInClust = {601, 1, 2, 0, 46};

  // set energy styles
  hEnePar   -> SetLineColor( get<0>(stylePar) );
  hEnePar   -> SetLineStyle( get<1>(stylePar) );
  hEnePar   -> SetLineWidth( get<2>(stylePar) );
  hEnePar   -> SetFillColor( get<0>(stylePar) );
  hEnePar   -> SetFillStyle( get<3>(stylePar) );
  hEnePar   -> SetMarkerColor( get<0>(stylePar) );
  hEnePar   -> SetMarkerStyle( get<4>(stylePar) );
  hEneHit   -> SetLineColor( get<0>(styleHit) );
  hEneHit   -> SetLineStyle( get<1>(styleHit) );
  hEneHit   -> SetLineWidth( get<2>(styleHit) );
  hEneHit   -> SetFillColor( get<0>(styleHit) );
  hEneHit   -> SetFillStyle( get<3>(styleHit) );
  hEneHit   -> SetMarkerColor( get<0>(styleHit) );
  hEneHit   -> SetMarkerStyle( get<4>(styleHit) );
  hEneClust -> SetLineColor( get<0>(styleClust) );
  hEneClust -> SetLineStyle( get<1>(styleClust) );
  hEneClust -> SetLineWidth( get<2>(styleClust) );
  hEneClust -> SetFillColor( get<0>(styleClust) );
  hEneClust -> SetFillStyle( get<3>(styleClust) );
  hEneClust -> SetMarkerColor( get<0>(styleClust) );
  hEneClust -> SetMarkerStyle( get<4>(styleClust) );
  hEneLead  -> SetLineColor( get<0>(styleLead) );
  hEneLead  -> SetLineStyle( get<1>(styleLead) );
  hEneLead  -> SetLineWidth( get<2>(styleLead) );
  hEneLead  -> SetFillColor( get<0>(styleLead) );
  hEneLead  -> SetFillStyle( get<3>(styleLead) );
  hEneLead  -> SetMarkerColor( get<0>(styleLead) );
  hEneLead  -> SetMarkerStyle( get<4>(styleLead) );
  hEneSum   -> SetLineColor( get<0>(styleSum) );
  hEneSum   -> SetLineStyle( get<1>(styleSum) );
  hEneSum   -> SetLineWidth( get<2>(styleSum) );
  hEneSum   -> SetFillColor( get<0>(styleSum) );
  hEneSum   -> SetFillStyle( get<3>(styleSum) );
  hEneSum   -> SetMarkerColor( get<0>(styleSum) );
  hEneSum   -> SetMarkerStyle( get<4>(styleSum) );

  // set difference styles
  hDiffClust -> SetLineColor( get<0>(styleClust) );
  hDiffClust -> SetLineStyle( get<1>(styleClust) );
  hDiffClust -> SetLineWidth( get<2>(styleClust) );
  hDiffClust -> SetFillColor( get<0>(styleClust) );
  hDiffClust -> SetFillStyle( get<3>(styleClust) );
  hDiffClust -> SetMarkerColor( get<0>(styleClust) );
  hDiffClust -> SetMarkerStyle( get<4>(styleClust) );
  hDiffLead  -> SetLineColor( get<0>(styleLead) );
  hDiffLead  -> SetLineStyle( get<1>(styleLead) );
  hDiffLead  -> SetLineWidth( get<2>(styleLead) );
  hDiffLead  -> SetFillColor( get<0>(styleLead) );
  hDiffLead  -> SetFillStyle( get<3>(styleLead) );
  hDiffLead  -> SetMarkerColor( get<0>(styleLead) );
  hDiffLead  -> SetMarkerStyle( get<4>(styleLead) );
  hDiffSum   -> SetLineColor( get<0>(styleSum) );
  hDiffSum   -> SetLineStyle( get<1>(styleSum) );
  hDiffSum   -> SetLineWidth( get<2>(styleSum) );
  hDiffSum   -> SetFillColor( get<0>(styleSum) );
  hDiffSum   -> SetFillStyle( get<3>(styleSum) );
  hDiffSum   -> SetMarkerColor( get<0>(styleSum) );
  hDiffSum   -> SetMarkerStyle( get<4>(styleSum) );

  // set multiplicity style
  hMultHit        -> SetLineColor( get<0>(styleHit) );
  hMultHit        -> SetLineStyle( get<1>(styleHit) );
  hMultHit        -> SetLineWidth( get<2>(styleHit) );
  hMultHit        -> SetFillColor( get<0>(styleHit) );
  hMultHit        -> SetFillStyle( get<3>(styleHit) );
  hMultHit        -> SetMarkerColor( get<0>(styleHit) );
  hMultHit        -> SetMarkerStyle( get<4>(styleHit) );
  hMultClust      -> SetLineColor( get<0>(styleClust) );
  hMultClust      -> SetLineStyle( get<1>(styleClust) );
  hMultClust      -> SetLineWidth( get<2>(styleClust) );
  hMultClust      -> SetFillColor( get<0>(styleClust) );
  hMultClust      -> SetFillStyle( get<3>(styleClust) );
  hMultClust      -> SetMarkerColor( get<0>(styleClust) );
  hMultClust      -> SetMarkerStyle( get<4>(styleClust) );
  hMultHitInClust -> SetLineColor( get<0>(styleInClust) );
  hMultHitInClust -> SetLineStyle( get<1>(styleInClust) );
  hMultHitInClust -> SetLineWidth( get<2>(styleInClust) );
  hMultHitInClust -> SetFillColor( get<0>(styleInClust) );
  hMultHitInClust -> SetFillStyle( get<3>(styleInClust) );
  hMultHitInClust -> SetMarkerColor( get<0>(styleInClust) );
  hMultHitInClust -> SetMarkerStyle( get<4>(styleInClust) );

  // legend parameters
  const pair<uint32_t, uint32_t> fLegCol = {0,   0};
  const pair<uint32_t, uint32_t> fLegSty = {0,   1};
  const pair<float,    float>    xLeg    = {0.1, 0.3};
  const pair<float,    float>    yLeg    = {0.1, 0.3};

  // create legends
  TLegend* lEne = new TLegend(xLeg.first, yLeg.first, xLeg.second, yLeg.second);
  lEne -> SetFillColor(fLegCol.first);
  lEne -> SetFillStyle(fLegSty.first);
  lEne -> SetLineColor(fLegCol.second);
  lEne -> SetLineStyle(fLegSty.second);
  lEne -> AddEntry(hEnePar.GetPtr(),   "Simulated Particle", "pf");
  lEne -> AddEntry(hEneHit.GetPtr(),   "BHCal Tiles", "pf");
  lEne -> AddEntry(hEneClust.GetPtr(), "BHCal Clusters", "pf");
  lEne -> AddEntry(hEneLead.GetPtr(),  "BHCal Lead Cluster", "pf");
  lEne -> AddEntry(hEneSum.GetPtr(),   "Sum of BHCal Clusters", "pf");

  TLegend* lDiff = new TLegend(xLeg.first, yLeg.first, xLeg.second, yLeg.second);
  lDiff -> SetFillColor(fLegCol.first);
  lDiff -> SetFillStyle(fLegSty.first);
  lDiff -> SetLineColor(fLegCol.second);
  lDiff -> SetLineStyle(fLegSty.second);
  lDiff -> AddEntry(hDiffClust.GetPtr(), "BHCal Clusters", "pf");
  lDiff -> AddEntry(hDiffLead.GetPtr(),  "BHCal Lead Cluster", "pf");
  lDiff -> AddEntry(hDiffSum.GetPtr(),   "Sum of BHCal Clusters", "pf");

  TLegend* lMult = new TLegend(xLeg.first, yLeg.first, xLeg.second, yLeg.second);
  lMult -> SetFillColor(fLegCol.first);
  lMult -> SetFillStyle(fLegSty.first);
  lMult -> SetLineColor(fLegCol.second);
  lMult -> SetLineStyle(fLegSty.second);
  lMult -> AddEntry(hMultHit.GetPtr(),        "BHCal Tiles", "pf");
  lMult -> AddEntry(hMultClust.GetPtr(),      "BHCal Clusters", "pf");
  lMult -> AddEntry(hMultHitInClust.GetPtr(), "BHCal Tiles per Cluster", "pf");

  // canvas parameters
  const uint32_t width  = 750;
  const uint32_t height = 750;
  const uint8_t  logY   = 1;
  const uint8_t  logZ   = 1;

  TCanvas *cEne = new TCanvas("cEne", "", width, height);
  cEne      -> cd();
  cEne      -> SetLogy(logY);
  hEnePar   -> Draw("hist");
  hEnePar   -> Draw("same");
  hEneHit   -> Draw("hist same");
  hEneHit   -> Draw("same");
  hEneClust -> Draw("hist same");
  hEneClust -> Draw("same");
  hEneLead  -> Draw("hist same");
  hEneLead  -> Draw("same");
  hEneSum   -> Draw("hist same");
  hEneSum   -> Draw("same");
  lEne      -> Draw();
  fOutput   -> cd();
  cEne      -> SaveAs("BHCalHitVsClusterEnergy.png");
  cEne      -> Write();
  cEne      -> Close();

  TCanvas *cDiff = new TCanvas("cDiff", "", width, height);
  cDiff      -> cd();
  cDiff      -> SetLogy(logY);
  hDiffClust -> Draw("hist");
  hDiffClust -> Draw("same");
  hDiffLead  -> Draw("hist same");
  hDiffLead  -> Draw("same");
  hDiffSum   -> Draw("hist same");
  hDiffSum   -> Draw("same");
  lDiff      -> Draw();
  fOutput    -> cd();
  cDiff      -> SaveAs("BHCalClustVsParDiff.png");
  cDiff      -> Write();
  cDiff      -> Close();

  TCanvas *cMult = new TCanvas("cMult", "", width, height);
  cMult           -> cd();
  cMult           -> SetLogy(logY);
  hMultHit        -> Draw("hist");
  hMultHit        -> Draw("same");
  hMultClust      -> Draw("hist same");
  hMultClust      -> Draw("same");
  hMultHitInClust -> Draw("same");
  hMultHitInClust -> Draw("hist same");
  lMult           -> Draw();
  fOutput         -> cd();
  cMult           -> SaveAs("BHCalTileVsClustMultiplicity.png");
  cMult           -> Write();
  cMult           -> Close();

  TCanvas *cSumVsFrac = new TCanvas("cSumVsFrac", "", width, height);
  cSumVsFrac -> cd();
  cSumVsFrac -> SetLogz(logZ);
  hSumVsFrac -> Draw("colz");
  fOutput    -> cd();
  cSumVsFrac -> SaveAs("BHCalVsBECalEnergy.png");
  cSumVsFrac -> Write();
  cSumVsFrac -> Close();

  // save output & exit -------------------------------------------------------

  // save output
  fOutput         -> cd();
  hEnePar         -> Write();
  hEneHit         -> Write();
  hEneClust       -> Write();
  hEneLead        -> Write();
  hEneSum         -> Write();
  hDiffClust      -> Write();
  hDiffLead       -> Write();
  hDiffSum        -> Write();
  hMultHit        -> Write();
  hMultClust      -> Write();
  hMultHitInClust -> Write();
  hSumVsFrac      -> Write();

  // close files and succesfully exit macro
  fOutput -> cd();
  fOutput -> Close();
  return 0;

}

// end ------------------------------------------------------------------------
