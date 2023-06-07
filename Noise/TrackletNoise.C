
#if !defined(__CLING__) || defined(__ROOTCLING__)


#endif

#include "DataFormatsTRD/Constants.h"

using namespace std::placeholders;

ROOT::RDF::RNode BuildMcmDF()
{

  // cout << "Creating RDataFrame" << endl;
  auto df = ROOT::RDataFrame(o2::trd::constants::NCHANNELSTOTAL / 21)
               .Define("sector", "rdfentry_ * 21 / o2::trd::constants::NCHANNELSPERSECTOR")
               .Define("layer", "(rdfentry_ * 21 % o2::trd::constants::NCHANNELSPERSECTOR) / o2::trd::constants::NCHANNELSPERLAYER")
               .Define("row", "(rdfentry_ * 21 % o2::trd::constants::NCHANNELSPERLAYER) / o2::trd::constants::NCHANNELSPERROW")
               .Define("col", "10 + rdfentry_ * 21 % o2::trd::constants::NCHANNELSPERROW")
               .Define("mcmcolglb", "(sector * o2::trd::constants::NCHANNELSPERROW + col) / o2::trd::constants::NADCMCM");

  return df;

}



// int getGlobalMCMColumn(sector)
// int side = tracklets[currenttracklet].getHCID() % 2; // 0: A-side, 1: B-side
// int stack = (detector % 30) / 6;
// int sec = detector / 30;
// int rowGlb = stack < 3 ? tracklets[currenttracklet].getPadRow() + stack * 16 : tracklets[currenttracklet].getPadRow() + 44 + (stack - 3) * 16; // pad row within whole sector
// int colGlb = tracklets[currenttracklet].getColumn() + sec * 8 + side * 4;
// mLayers[layer]->Fill(rowGlb, colGlb);

struct mcminfo_t {
  int detector{-1};
  int rob{-1};
  int mcm{-1};

  uint32_t noisychannels{0};
  uint32_t adcmask{0};
  size_t ntracklets{0};
};

class TrackletMCMStats
{
public:
  TrackletMCMStats(const long timestamp, const char* url = 0);
  float GetNTrackletsPerMCMForChannel(ULong64_t sector, ULong64_t layer, ULong64_t row, ULong64_t col);
  float GetNTrackletsPerMCM(ULong64_t mcmindex);

  /// Find the largest number of tracklets for any MCM
  int GetMaxNTracklets();

  ROOT::RDF::RNode AddToDF(ROOT::RDF::RNode df);
  map<int, mcminfo_t> FindNoisyMCMs(int min_tracklets);

  // protected:
  TH2* mQCLayerHistos[6];
};

float TrackletMCMStats::GetNTrackletsPerMCM(ULong64_t mcmindex)
{
  int col = mcmindex % 8;
  mcmindex /= 8;

  int row = mcmindex % 76;
  mcmindex /= 76;

  int layer = mcmindex % 6;
  int sector = mcmindex / 6;

  int colglb = sector*8 + col;

  return mQCLayerHistos[layer]->GetBinContent(row + 1, colglb + 1);
}

float TrackletMCMStats::GetNTrackletsPerMCMForChannel(ULong64_t sector, ULong64_t layer, ULong64_t row, ULong64_t col)
{
  int colglb = (sector * o2::trd::constants::NCHANNELSPERROW + col) / o2::trd::constants::NADCMCM;
  return mQCLayerHistos[layer]->GetBinContent(row + 1, colglb + 1);
}

int TrackletMCMStats::GetMaxNTracklets()
{
  int max = 0;
  for (int l=0; l<6; ++l) {
    for (int r=1; r<=76; ++r) {
      for (int c=1; c<=18*8; ++c) {
        if (mQCLayerHistos[l]->GetBinContent(r, c) > max) {
          max = mQCLayerHistos[l]->GetBinContent(r, c);
        }
      }
    }
  }
  return max;
}

ROOT::RDF::RNode TrackletMCMStats::AddToDF(ROOT::RDF::RNode df)
{
  return df
    .Define("ntracklets", [this](ULong64_t sector, ULong64_t layer, ULong64_t row, ULong64_t col) {
      return this->GetNTrackletsPerMCMForChannel(sector, layer, row, col);
    }, {"sector", "layer", "row", "col"});

  // P.S.: the above solution looks clumsy to me, and I hoped I could pass a pointer to the member
  // function, bound to this, to `Define`. I tried the following, but cave a compiler warning:
  // return df.Define("ntracklets", this->GetNTrackletsPerMCMForChannel, {"sector", "layer", "row", "col"});
}



map<int, mcminfo_t> TrackletMCMStats::FindNoisyMCMs(int min_tracklets)
{
  map<int, mcminfo_t> mcminfo;

  for(int sector=0; sector<18; ++sector) {
    for(int layer=0; layer<6; ++layer) {
      for(int row=0; row<76; ++row) {
        for (int col = 0; col < 168; col += 21) {
          int ntrkl = GetNTrackletsPerMCMForChannel(sector, layer, row, col);
          if (ntrkl>min_tracklets) {
            int mcmidx = ((sector*6 + layer)*76 + row)*8 + col/21;

            int det, rob, mcm, ch;
            o2::trd::HelperMethods::getPositionFromGlobalChannelIndex(mcmidx*21, det, rob, mcm, ch);

            mcminfo[mcmidx].detector = det;
            mcminfo[mcmidx].rob = rob;
            mcminfo[mcmidx].mcm = mcm;
            mcminfo[mcmidx].ntracklets = ntrkl;
          }
        }
      }
    }    
  }
  return mcminfo;
}

TrackletMCMStats::TrackletMCMStats(const long timestamp, const char* url)
{
  auto& ccdbmgr = o2::ccdb::BasicCCDBManager::instance();
  if (url != 0) {
    ccdbmgr.setURL(url);
  }
  ccdbmgr.setTimestamp(timestamp);

  for (int L=0; L<6; ++L) {
    mQCLayerHistos[L] = ccdbmgr.get<TH2F>(Form("qc/TRD/MO/TrackletsTask/TrackletsPerLayer/layer%d", L));
    if (!mQCLayerHistos[L]) {
      LOG(fatal) << "No ";
    }
  }
}



void TrackletNoise()
{
  // auto df = BuildMcmDF();

  // Initialize statistics object with tracklets per MCM from test CCDB, c.f.
  // http:// ccdb-test.cern.ch:8080/browse/qc/TRD/MO/TrackletsTask/TrackletsPerLayer/layer0
  TrackletMCMStats mcmstats(1684761589464, "http://ccdb-test.cern.ch:8080");

  // Draw the layer plots
  for (int L = 0; L < 6; ++L) {
    TCanvas* cnv = new TCanvas(Form("layer%d", L), Form("layer%d", L));
    cnv->SetLogz();
    mcmstats.mQCLayerHistos[L]->Draw("colz");
  }

  // Build a RDataFrame with one entry per MCM
  auto df = mcmstats.AddToDF(BuildMcmDF());

  // Draw the number of tracklets per MCM
  TCanvas* cnv = new TCanvas("ntracklets", "Tracklets per layer");
  cnv->SetLogy();
  auto h = df.Histo1D("ntracklets");
  h->DrawClone();

  for (auto& [idx, mcm] : mcmstats.FindNoisyMCMs(100e3)) {
    cout << mcm.detector << " " << mcm.rob << ":" << mcm.mcm << endl;
  }
}