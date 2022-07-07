
#include "../macros/HelperFunctions.C"
#include "../macros/FmtFormat.C"

#include <DataFormatsTRD/Digit.h>
#include <DataFormatsTRD/Tracklet64.h>
#include <DataFormatsTRD/TriggerRecord.h>
#include <DataFormatsTRD/Hit.h>
#include <DataFormatsTRD/Constants.h>

#include "ReconstructionDataFormats/TrackTPCITS.h"
#include <DataFormatsTPC/TrackTPC.h>
#include <CommonDataFormat/TFIDInfo.h>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TLine.h>

template <typename value_t, template <typename T> typename container_t>
struct myrange
{
  // typedef value_t value_type;
  typedef typename container_t<value_t>::iterator iterator;

  iterator &begin() { return b; }
  iterator &end() { return e; }
  iterator b, e;

  size_t length() { return e - b; }
  void sort(bool (*comp)(const value_t &a, const value_t &b))
  {
    std::stable_sort(b, e, comp);
  }
};

struct RawDataSpan
{
  myrange<o2::trd::Digit, TTreeReaderArray> digits;
  myrange<o2::trd::Tracklet64, TTreeReaderArray> tracklets;

  vector<o2::tpc::TrackTPC> tpctracks;

  pair<int, int> getMaxADCsumAndChannel();
  int getMaxADCsum(){ return getMaxADCsumAndChannel().first; }
};

pair<int, int> RawDataSpan::getMaxADCsumAndChannel()
{
  int maxch = -1, maxval = -1;
  for (auto digit : digits) {
    if (digit.getADCsum() > maxval) {
      maxval = digit.getADCsum();
      maxch = digit.getChannel();
    }
  // cout << digit << endl;
  }

  return make_pair(maxval, maxch);
}

ostream& operator<<(ostream& os, RawDataSpan& sp)
{
  os << "raw data span with " << sp.digits.length() << " digits and " 
     << sp.tracklets.length() << " tracklets";
  return os;
}

bool order_digit(const o2::trd::Digit &a, const o2::trd::Digit &b)
{
  if (a.getDetector() != b.getDetector())
    return a.getDetector() < b.getDetector();

  if (a.getPadRow() != b.getPadRow())
    return a.getPadRow() < b.getPadRow();

  if (a.getROB() != b.getROB())
    return a.getROB() < b.getROB();

  if (a.getMCM() != b.getMCM())
    return a.getMCM() < b.getMCM();

  if (a.getChannel() != b.getChannel())
    return a.getChannel() < b.getChannel();

  return true;
}

bool order_tracklet(const o2::trd::Tracklet64 &a, const o2::trd::Tracklet64 &b)
{
  // upper bits of hcid and padrow from Tracklet64 word
  const uint64_t det_row_mask = 0x0ffde00000000000;

  // lowest bit of hcid (side), MCM col and pos from Tracklet64 word
  const uint64_t col_pos_mask = 0x00011fff00000000;

  auto a_det_row = a.getTrackletWord() & det_row_mask;
  auto b_det_row = b.getTrackletWord() & det_row_mask;

  if (a_det_row != b_det_row)
    return a_det_row < b_det_row;

  auto a_col_pos = a.getTrackletWord() & col_pos_mask;
  auto b_col_pos = b.getTrackletWord() & col_pos_mask;

  return a_col_pos < b_col_pos;  
};


struct ClassifierByDetector
{
  static uint32_t key(const o2::trd::Digit &x) { return x.getDetector(); }
  static uint32_t key(const o2::trd::Tracklet64 &x) { return x.getDetector(); }

  static bool comp_digits(const o2::trd::Digit &a, const o2::trd::Digit &b)
  {
    return key(a) != key(b);
  }

  static bool comp_tracklets(const o2::trd::Tracklet64 &a, const o2::trd::Tracklet64 &b)
  {
    return key(a) != key(b);
  }
};

struct ClassifierByPadRow
{
  static uint32_t key(const o2::trd::Digit &x) { 
    return 100*x.getDetector() + x.getPadRow();
  }
  static uint32_t key(const o2::trd::Tracklet64 &x) {
    return 100 * x.getDetector() + x.getPadRow();
  }

  static bool comp_digits(const o2::trd::Digit &a, const o2::trd::Digit &b)
  {
    return key(a) != key(b);
  }

  static bool comp_tracklets(const o2::trd::Tracklet64 &a, const o2::trd::Tracklet64 &b)
  {
    return key(a) != key(b);
  }
};

struct ClassifierByMCM
{
  template<typename T>
  static uint32_t key(const T &x) { 
    return 1000*x.getDetector() + 10*x.getPadRow() + 4*(x.getROB()%2) + x.getMCM()%4;
  }

  static bool comp_digits(const o2::trd::Digit &a, const o2::trd::Digit &b)
  {
    return key(a) != key(b);
  }

  static bool comp_tracklets(const o2::trd::Tracklet64 &a, const o2::trd::Tracklet64 &b)
  {
    return key(a) != key(b);
  }
};

// struct ClassifierByContinousRegion
// {
//   template<typename T>
//   static uint32_t mcmid(const T &x) { 
//     return 1000*x.getDetector() + 10*x.getPadRow() + 4*(x.getROB()%2) + x.getMCM()%4;
//   }

//   template<typename T>
//   static uint32_t mcmid(const T &x) { 
//     return 1000*x.getDetector() + 10*x.getPadRow() + 4*(x.getROB()%2) + x.getMCM()%4;
//   }

//   static bool comp_digits(const o2::trd::Digit &a, const o2::trd::Digit &b)
//   {
//     if (mcmid(a) != mcmid(b)) {
//       return false;
//     } else if ( abs(a.getChannel()-b.getChannel()) == 1 ) {
//       return true;
//     } else {
//       return false;
//     }
//   }

//   // static bool comp_tracklets(const o2::trd::Tracklet64 &a, const o2::trd::Tracklet64 &b)
//   // {
//   //   return key(a) != key(b);
//   // }
// };


template<typename classifier>
class RawDataPartitioner : public map<uint32_t, RawDataSpan>
{
public:

  RawDataPartitioner(RawDataSpan event)
  {
    // sort digits and tracklets
    event.digits.sort(order_digit);
    event.tracklets.sort(order_tracklet);

    // add all the digits to a map that contains all the 
    for (auto cur = event.digits.b; cur != event.digits.e; /* noop */ ) {
      // let's look for pairs where the comparision function is true, i.e.
      // where two adjacent elements do not have equal keys
      auto nxt = std::adjacent_find(cur, event.digits.e, classifier::comp_digits);
      // if adjacent_find found a match, it returns the first element of the
      // adjacent pair, but we need the second one -> increment the result
      if (nxt != event.digits.e)
      {
        ++nxt;
      }
      // finally, we add the found elements to the map with the current range
      (*this)[classifier::key(*cur)].digits.b = cur; 
      (*this)[classifier::key(*cur)].digits.e = nxt;
      cur = nxt;
    }

    for (auto cur = event.tracklets.b; cur != event.tracklets.e; /* noop */) {
      auto nxt = std::adjacent_find(cur, event.tracklets.e, classifier::comp_tracklets);
      if (nxt != event.tracklets.e) {
        ++nxt;
      }
      (*this)[classifier::key(*cur)].tracklets.b = cur; 
      (*this)[classifier::key(*cur)].tracklets.e = nxt;
      cur = nxt;
    }
    cout << "Found " << this->size() << " spans" << endl;
  }
};

// ========================================================================
// the DataManager class
// ========================================================================


class DataManager
{

public:
  DataManager(std::string dir = "data/");

  void SetMatchWindowTPC(float min, float max)
  { mMatchTimeMinTPC=min; mMatchTimeMaxTPC=max; }

  bool NextTimeFrame();
  bool NextEvent();

  // access time frame info
  o2::dataformats::TFIDInfo GetTimeFrameInfo();
  TTreeReaderArray<o2::tpc::TrackTPC> *GetTimeFrameTPCTracks()
  {return tpctracks; }
  
  TTreeReaderArray<o2::dataformats::TrackTPCITS> *GetTimeFrameTracks()
  { return tracks; }

  // access event info
  RawDataSpan GetEvent();
  float GetTriggerTime();

  size_t GetTimeFrameNumber() { return iTimeframe; }
  size_t GetEventNumber() { return iEvent; }

private: 
  TFile *mainfile;
  TTree* datatree; // tree and friends from digits, tracklets files
  TTreeReader *datareader;

  TTreeReaderArray<o2::trd::Hit>* hits;
  TTreeReaderArray<o2::trd::Digit>* digits;
  TTreeReaderArray<o2::trd::Tracklet64>* tracklets;
  TTreeReaderArray<o2::trd::TriggerRecord>* trgrecords;

  TTreeReaderArray<o2::dataformats::TrackTPCITS> *tracks;
  TTreeReaderArray<o2::tpc::TrackTPC> *tpctracks;

  o2::trd::TriggerRecord mTriggerRecord;

  std::vector<o2::dataformats::TFIDInfo> *tfids;

  size_t iTimeframe, iEvent;
  float mMatchTimeMinTPC{-10.0}, mMatchTimeMaxTPC{20.0};

  template <typename T>
  TTreeReaderArray<T> *AddReaderArray(std::string file, std::string tree, std::string branch);

  template <typename T>
  void AddReaderArray(TTreeReaderArray<T> *& array, std::string file, std::string tree, std::string branch);
};

DataManager::DataManager(std::string dir)
    : mainfile(0), datatree(0), datareader(0),
      hits(0), digits(0), tracklets(0), trgrecords(0), tpctracks(0),
      tfids(0),
      iTimeframe(0), iEvent(0)
{

  // We allways need the trigger records, which are stored in trdtracklets.root.
  // While at it, let's also set up reading the tracklets.
  mainfile = new TFile((dir+"trdtracklets.root").c_str());
  mainfile->GetObject("o2sim", datatree);
  datareader = new TTreeReader(datatree);

  // set up the branches we want to read
  tracklets = new TTreeReaderArray<o2::trd::Tracklet64>(*datareader, "Tracklet");
  trgrecords = new TTreeReaderArray<o2::trd::TriggerRecord>(*datareader, "TrackTrg");

  // For data, we need info about time frames to match ITS and TPC tracks to 
  // trigger records.
  TFile* fInTFID = TFile::Open((dir+"o2_tfidinfo.root").c_str());
  if (fInTFID) {
    // for a simulation this file is not available
    tfids = (std::vector<o2::dataformats::TFIDInfo> *)fInTFID->Get("tfidinfo");
  }

  // Let's try to add other data
  // digits = AddReaderArray<o2::trd::Digit>(dir+"trddigits.root", "o2sim", "TRDDigit");
  // tpctracks = AddReaderArray<o2::tpc::TrackTPC>(dir + "tpctracks.root", "tpcrec", "TPCTracks");
  AddReaderArray(digits, dir+"trddigits.root", "o2sim", "TRDDigit");
  AddReaderArray(tpctracks, dir + "tpctracks.root", "tpcrec", "TPCTracks");
  AddReaderArray(tracks, dir + "o2match_itstpc.root", "matchTPCITS", "TPCITS");

  // ConnectMCHitsFile(dir+"o2sim_HitsTRD.root");
}

template <typename T>
TTreeReaderArray<T> *DataManager::AddReaderArray(std::string file, std::string tree, std::string branch)
{
  if (gSystem->AccessPathName(file.c_str())) {
    // file was not found
    return nullptr;
  }

  // the file exists, everything else should work
  datatree->AddFriend(tree.c_str(), file.c_str());
  return new TTreeReaderArray<T>(*datareader, branch.c_str());
}

template <typename T>
void DataManager::AddReaderArray(TTreeReaderArray<T> *&array, std::string file, std::string tree, std::string branch)
{
  if (gSystem->AccessPathName(file.c_str())) {
    // file was not found
    return;
  }

  // the file exists, everything else should work
  datatree->AddFriend(tree.c_str(), file.c_str());
  array = new TTreeReaderArray<T>(*datareader, branch.c_str());
}

bool DataManager::NextTimeFrame()
{
  if (datareader->Next())
  {
    iEvent = 0;
    iTimeframe++;
    cout << "## Time frame " << iTimeframe << endl;

    for (auto &tracklet : *tracklets) {
      tracklet.setPosition(tracklet.getPosition() ^ 0x80);
      tracklet.setSlope(tracklet.getSlope() ^ 0x80);
    }

    return true;
  }
  else
  {
    return false;
  }
}

bool DataManager::NextEvent()
{
  // get the next trigger record
  if (iEvent >= trgrecords->GetSize()) {
    return false;
  }

  // load the hits for the next event
  // if ( ! rdrhits->Next() ) {
  //   cout << "no hits found for event " << tfno << ":" << iEvent << endl;
  //   return false;
  // }

  mTriggerRecord = trgrecords->At(iEvent);
  // cout << mTriggerRecord << endl;

  cout << "## Event " << iTimeframe << ":" << iEvent << ":  "
       //  << hits->GetSize() << " hits   "
       << mTriggerRecord.getNumberOfDigits() << " digits and "
       << mTriggerRecord.getNumberOfTracklets() << " tracklets" << endl;

  iEvent++;
  return true;
}

RawDataSpan DataManager::GetEvent()
{
  RawDataSpan ev;
  // ev.digits = GetDigits();
  ev.digits.b = digits->begin() + mTriggerRecord.getFirstDigit();
  ev.digits.e = ev.digits.begin() + mTriggerRecord.getNumberOfDigits();
  ev.tracklets.b = tracklets->begin() + mTriggerRecord.getFirstTracklet();
  ev.tracklets.e = ev.tracklets.begin() + mTriggerRecord.getNumberOfTracklets();

  auto evtime = GetTriggerTime();
  if (tpctracks) {
    for (auto &track : *tpctracks) {
      //   // auto tracktime = track.getTimeMUS().getTimeStamp();
      auto dtime = track.getTime0() / 5.0 - evtime;
      if (dtime > mMatchTimeMinTPC && dtime < mMatchTimeMaxTPC) {
        ev.tpctracks.push_back(track);
      }
    }
  }

  return ev;
}

o2::dataformats::TFIDInfo DataManager::GetTimeFrameInfo()
{
  if (tfids) {
    return tfids->at(iTimeframe-1);
  } else {
    return o2::dataformats::TFIDInfo();
  }
}

float DataManager::GetTriggerTime()
{
  auto tfid = GetTimeFrameInfo();

  if (tfid.isDummy()) {
    return mTriggerRecord.getBCData().bc2ns() * 1e-3;
  } else {
    o2::InteractionRecord intrec = {0, tfid.firstTForbit};
    return mTriggerRecord.getBCData().differenceInBCMS(intrec);
  }
}
// void DataManager::ConnectMCHitsFile(std::string fname)
// {
//   // ----------------------------------------------------------------------
//   // set up data structures for reading

//   if (fhits || trhits) {
//     cerr << "Hits file seems to be connected." << endl;
//     return;
//   }

//   fhits = new TFile(fname.c_str());
//   fhits->GetObject("o2sim", trhits);

//   rdrhits = new TTreeReader(trhits);
//   hits = new TTreeReaderArray<o2::trd::Hit>(*rdrhits, "TRDHit");
// }

// ========================================================================
// Drawing routines
// ========================================================================

TVirtualPad *DrawPadRow(RawDataSpan &padrow, TVirtualPad *pad = NULL, TH2F *adcmap = NULL)
{
  auto x = *padrow.digits.begin();
  string desc = fmt::format("{:m}", x);
  string name = fmt::format("det{:03d}_rob{:d}_mcm{:02d}",
                            x.getDetector(), x.getROB(), x.getMCM());

  if (pad == NULL)
  {
    pad = new TCanvas(desc.c_str(), desc.c_str(), 1200, 500);
    pad->cd();
  }
  // else
  // {
  //   pad->SetName(name.c_str());
  //   pad->SetTitle(desc.c_str());
  // }

  if (adcmap == NULL) {
    adcmap = new TH2F(name.c_str(), (desc + ";pad;time bin").c_str(), 144, 0., 144., 30, 0., 30.);
  }

  for (auto digit : padrow.digits)
  {
    if (digit.isSharedDigit())
    {
      continue;
    }

    auto adc = digit.getADC();
    for (int tb = 0; tb < 30; ++tb)
    {
      adcmap->Fill(digit.getPadCol(), tb, adc[tb]);
    }
  }
  adcmap->SetStats(0);
  adcmap->Draw("colz");
  adcmap->Draw("text,same");

  TLine trkl;
  trkl.SetLineWidth(2);
  trkl.SetLineColor(kRed);
  // trkl.SetLineStyle(kDotted);

  TLine trkl2;
  trkl2.SetLineWidth(4);
  trkl2.SetLineStyle(kDashed);
  trkl2.SetLineColor(kBlack);

  for (auto tracklet : padrow.tracklets)
  {
    auto pos = PadPosition(tracklet);
    auto ypos = UncalibratedPad(tracklet);
    auto slope = Slope(tracklet);
    trkl.DrawLine(pos, 0, pos - 30*slope, 30);
    // trkl2.DrawLine(ypos, 0, ypos - 30*slope, 30);
  }

  return pad;
}

TPad *DrawMCM(RawDataSpan &mcm, TPad *pad)
{
  auto x = *mcm.digits.begin();
  string desc = fmt::format("{:m}", x);
  string name = fmt::format("det{:03d}_rob{:d}_mcm{:02d}",
                            x.getDetector(), x.getROB(), x.getMCM());

  if (pad == NULL)
  {
    pad = new TCanvas(desc.c_str(), desc.c_str(), 800, 600);
  }
  else
  {
    pad->SetName(name.c_str());
    pad->SetTitle(desc.c_str());
  }
  pad->cd();

  TH2F *digit_disp = new TH2F(desc.c_str(), (desc + ";ADC channel;time bin").c_str(), 21, 0., 21., 30, 0., 30.);

  for (auto digit : mcm.digits)
  {
    auto adc = digit.getADC();
    for (int tb = 0; tb < 30; ++tb)
    {
      digit_disp->Fill(digit.getChannel(), tb, adc[tb]);
    }
  }
  digit_disp->SetStats(0);
  digit_disp->Draw("colz");

  TLine trkl;
  trkl.SetLineColor(kRed);
  trkl.SetLineWidth(3);

  for (auto tracklet : mcm.tracklets)
  {
    auto pos = PadPositionMCM(tracklet);
    auto slope = Slope(tracklet);
    trkl.DrawLine(pos, 0, pos + 30 * slope, 30);
  }

  return pad;
}
