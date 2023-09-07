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

using namespace std;
using namespace ROOT;



// bhcal benchmarks -----------------------------------------------------------

int bhcal(std::string file) {

  // histogram binning
  const tuple<int, double, double> eneBins  = make_tuple(200, 0.,   100.);
  const tuple<int, double, double> diffBins = make_tuple(200, -10., 10.);
  const tuple<int, double, double> fracBins = make_tuple(100, 0.,   10.);
  const tuple<int, double, double> multBins = make_tuple(100, 0.,   100.);

  // histogram definitions
  const vector<RDF::TH1DModel> vecHistDefs1D = {
    RDF::TH1DModel("hEneHit",         "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    RDF::TH1DModel("hEneClust",       "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    RDF::TH1DModel("hEneLead",        "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    RDF::TH1DModel("hEneSum",         "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    RDF::TH1DModel("hDiffClust",      "", get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    RDF::TH1DModel("hDiffLead",       "", get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    RDF::TH1DModel("hDiffSum",        "", get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    RDF::TH1DModel("hFracBHxECal",    "", get<0>(fracBins), get<1>(fracBins), get<2>(fracBins)),
    RDF::TH1DModel("hMultHit",        "", get<0>(multBins), get<1>(multBins), get<2>(multBins)),
    RDF::TH1DModel("hMultClust",      "", get<0>(multBins), get<1>(multBins), get<2>(multBins)),
    RDF::TH1DModel("hMultHitInClust", "", get<0>(multBins), get<2>(multBins), get<2>(multBins))
  };
  const vector<RDF::TH2DModel> vecHistDefs2D = {
    RDF::TH2DModel("hEneDiffClustVsParticle", "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins),  get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    RDF::TH2DModel("hEneDiffLeadVsParticle",  "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins),  get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    RDF::TH2DModel("hEneDiffSumVsParticle",   "", get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins),  get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    RDF::TH2DModel("hEneDiffVsFracBHxECal",   "", get<0>(fracBins), get<1>(fracBins), get<2>(fracBins), get<0>(diffBins), get<1>(diffBins), get<1>(diffBins))
  };

  // enable multithreading
  EnableImplicitMT(kNumThreads);

  // open dataframe
  RDataFrame frame("events", file);

  // make sure file isn't empty
  const uint64_t events = frame.Count();
  if (events == 0) {
    cerr << "Error: No events found!" << endl;
    assert(events > 0);
  }

  /* TODO operations go here */

  // succesfully exit macro
  return 0;

}

// end ------------------------------------------------------------------------
