// ----------------------------------------------------------------------------
// 'bhcal.cxx'
// Derek Anderson
// 09.07.2023
//
// ePIC BHCal benchmark macro.
// ----------------------------------------------------------------------------

// c utilities
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>
// root classes
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
// dataframe related classes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
// benchmark utilities
#include "common_bench/mt.h"
#include "common_bench/util.h"
#include "common_bench/plot.h"
#include "common_bench/benchmark.h"
// formatting utilities
#include "fmt/core.h"
#include "fmt/color.h"
// misc
#include "nlohmann/json.hpp"

// make common namespaces implicit
using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps; 

// set up aliases
using TH1Def = ROOT::RDF::TH1DModel;
using TH2Def = ROOT::RDF::TH2DModel;



// bhcal benchmarks -----------------------------------------------------------

int bhcal(std::string file) {

  // enable multithreading
  EnableImplicitMT(kNumThreads);

  // open dataframe
  RDataFrame frame("events", file);

  // make sure file isn't empty
  auto events = frame.Count();
  if (events == 0) {
    cerr << "Error: No events found!" << endl;
    assert(events > 0);
  }

  // lambdas
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

  // TODO make into template
  auto getMultiplicity = [](const RVec<float> &collection) {
    return collection.size();
  };

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

  // histogram binning
  const tuple<int, double, double> eneBins  = make_tuple(200, 0.,   100.);
  const tuple<int, double, double> diffBins = make_tuple(200, -10., 10.);
  const tuple<int, double, double> fracBins = make_tuple(100, 0.,   10.);
  const tuple<int, double, double> multBins = make_tuple(100, 0.,   100.);

  // histogram definitions
  const vector<TH1Def> vecHistDefs1D = {
    TH1Def("hEneHit",         "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneClust",       "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneLead",        "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneSum",         "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hDiffClust",      "", get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hDiffLead",       "", get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hDiffSum",        "", get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hFracBHxECal",    "", get<0>(fracBins), get<1>(fracBins), get<2>(fracBins)),
    TH1Def("hMultHit",        "", get<0>(multBins), get<1>(multBins), get<2>(multBins)),
    TH1Def("hMultClust",      "", get<0>(multBins), get<1>(multBins), get<2>(multBins)),
    TH1Def("hMultHitInClust", "", get<0>(multBins), get<2>(multBins), get<2>(multBins))
  };
  const vector<TH2Def> vecHistDefs2D = {
    TH2Def("hEneDiffClustVsParticle", "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins),  get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH2Def("hEneDiffLeadVsParticle",  "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins),  get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH2Def("hEneDiffSumVsParticle",   "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins),  get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH2Def("hEneDiffVsFracBHxECal",   "", get<0>(fracBins), get<1>(fracBins), get<2>(fracBins), get<0>(diffBins), get<1>(diffBins), get<1>(diffBins))
  };

  // get 1D histograms
  // TODO collect into vector
  auto hEnePar         = analysis.Histo1D(vecHistDefs1D[0],  "enePar");
  auto hEneHit         = analysis.Histo1D(vecHistDefs1D[1],  "HcalBarrelRecHits.energy");
  auto hEneClust       = analysis.Histo1D(vecHistDefs1D[2],  "HcalBarrelClusters.energy");
  auto hEneLead        = analysis.Histo1D(vecHistDefs1D[3],  "eneLeadHCal");
  auto hEneSum         = analysis.Histo1D(vecHistDefs1D[4],  "eneSumHCal");
  auto hDiffClust      = analysis.Histo1D(vecHistDefs1D[5],  "diffClustHCal");
  auto hDiffLead       = analysis.Histo1D(vecHistDefs1D[6],  "diffLeadHCal");
  auto hDiffSum        = analysis.Histo1D(vecHistDefs1D[7],  "diffSumHCal");
  auto hFracHxECal     = analysis.Histo1D(vecHistDefs1D[8],  "fracSumHxECal");
  auto hMultHit        = analysis.Histo1D(vecHistDefs1D[9],  "multHit");
  auto hMultClust      = analysis.Histo1D(vecHistDefs1D[10], "multClust");
  auto hMultHitInClust = analysis.Histo1D(vecHistDefs1D[11], "HcalBarrelClusters.nhits");

  // get 2D histograms
  // TODO collect into vector
  auto hDiffClustVsParticle = analysis.Histo2D(vecHistDefs2D[0], "enePar",        "diffClustHCal");
  auto hDiffLeadVsParticle  = analysis.Histo2D(vecHistDefs2D[1], "enePar",        "diffLeadHCal");
  auto hDiffSumVsParticle   = analysis.Histo2D(vecHistDefs2D[2], "enePar",        "diffSumHCal");
  auto hSumVsFracHxECal     = analysis.Histo2D(vecHistDefs2D[3], "fracSumHxECal", "eneSumHCal");

  /* TODO make plots */

  // succesfully exit macro
  return 0;

}

// end ------------------------------------------------------------------------
