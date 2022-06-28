
#include "../macros/HelperFunctions.C"
#include "../macros/FmtFormat.C"

#include <DataFormatsTRD/Digit.h>
#include <DataFormatsTRD/Tracklet64.h>
#include <DataFormatsTRD/TriggerRecord.h>
#include <DataFormatsTRD/Hit.h>
#include <DataFormatsTRD/Constants.h>

#include <DataFormatsTPC/TrackTPC.h>

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
};

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


class DataManager
{

public:

  DataManager(std::string dir="data/")
  : mainfile(0), datatree(0), datareader(0),
    hits(0), digits(0), tracklets(0), trgrecords(0),
    iTimeframe(0), iEvent(0)
  {

    std::string testfiles[] = { "trddigits.root", "trdtracklets.root" };

    for ( auto fname : testfiles ) {
      auto fullname = dir + fname;
      if(gSystem->AccessPathName(fullname.c_str())) { 
        cout << "WARN: The file " << fname << " doesn't exist. Skipping." 
             << endl;
        continue;
      } 

      if ( ! datatree ) {
        mainfile = new TFile(fullname.c_str());
        mainfile->GetObject("o2sim", datatree);
        datareader = new TTreeReader(datatree);
      } else {
        datatree->AddFriend("o2sim", fullname.c_str());
      }
    }

    // auto tracksfilename = dir + "tpctracks.root";
    // if (gSystem->AccessPathName(tracksfilename.c_str())) {
    //   datatree->AddFriend("tpcrec", tracksfilename.c_str());
    // }
     
    // set up the branches we want to read
    digits = new TTreeReaderArray<o2::trd::Digit>(*datareader, "TRDDigit");
    tracklets = new TTreeReaderArray<o2::trd::Tracklet64>(*datareader, "Tracklet");
    trgrecords = new TTreeReaderArray<o2::trd::TriggerRecord>(*datareader, "TrackTrg");
    // tpctracks = new TTreeReaderArray<o2::tpc::TrackTPC>(*datareader, "TPCTracks");

    // ConnectMCHitsFile(dir+"o2sim_HitsTRD.root");
    }

  bool NextTimeFrame()
  {
    if (datareader->Next() ) {
      iEvent = 0;
      iTimeframe++;
      cout << "## Time frame " << iTimeframe << endl;

      for (auto &tracklet : *tracklets)
      {
        tracklet.setPosition(tracklet.getPosition() ^ 0x80);
        tracklet.setSlope(tracklet.getSlope() ^ 0x80);
      }

      // cout << "Found " << tpctracks->GetSize() << " TPC tracks" << endl;

      return true;
    } else {
      return false;
    }
  }

  bool NextEvent()
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

  RawDataSpan GetEvent()
  {
    RawDataSpan ev;
    // ev.digits = GetDigits();
    ev.digits.b = digits->begin() + mTriggerRecord.getFirstDigit();
    ev.digits.e = ev.digits.begin() + mTriggerRecord.getNumberOfDigits();
    ev.tracklets.b = tracklets->begin() + mTriggerRecord.getFirstTracklet();
    ev.tracklets.e = ev.tracklets.begin() + mTriggerRecord.getNumberOfTracklets();
    return ev;
  }

  size_t GetTimeFrameNumber() { return iTimeframe; }
  size_t GetEventNumber() { return iEvent; }


protected :
      // void ConnectMCHitsFile(std::string fname)
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

private: 
  TFile *mainfile;
  // TFile *fhits, *fdigits, *ftracklets;
  // TTree *trhits, *trdigits, *trtracklets, *trtrgrec;
  TTree* datatree; // tree and friends from digits, tracklets files
  // TTreeReader *rdrhits, *rdrdigits, *rdrtracklets, *rdrtrgreg;
  TTreeReader *datareader;

  TTreeReaderArray<o2::trd::Hit>* hits;
  TTreeReaderArray<o2::trd::Digit>* digits;
  TTreeReaderArray<o2::trd::Tracklet64>* tracklets;
  TTreeReaderArray<o2::trd::TriggerRecord>* trgrecords;

  TTreeReaderArray<o2::tpc::TrackTPC> *tpctracks;

  // DigitRange ev_digits;//, sel_digits;
  o2::trd::TriggerRecord mTriggerRecord;

  size_t iTimeframe, iEvent;
};

TVirtualPad *DrawPadRow(RawDataSpan &padrow, TVirtualPad *pad = NULL, TH2F* adcmap = NULL)
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
