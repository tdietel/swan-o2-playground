
#include "../macros/HelperFunctions.C"
#include "../macros/FmtFormat.C"

#include <DataFormatsTRD/Digit.h>
#include <DataFormatsTRD/Tracklet64.h>
#include <DataFormatsTRD/TriggerRecord.h>
#include <DataFormatsTRD/Hit.h>
#include <DataFormatsTRD/Constants.h>

#include <TRDBase/Geometry.h>

#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <ReconstructionDataFormats/TrackTPCITS.h>
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
#include <TMarker.h>

#include <vector>

using namespace std;

// ========================================================================
// ========================================================================
//
// ChamberSpacePoint:  A position in spatial (x,y,z) and raw/digit
// coordinates (det,row,col,tb).
//
// ========================================================================
// ========================================================================

class ChamberSpacePoint
{
 public:
  ChamberSpacePoint(int det = -999) : mDetector(det){};
  ChamberSpacePoint(o2::track::TrackParCov& t, o2::trd::Geometry* geo = nullptr);

  bool isValid() const { return mDetector >= 0; }
  float getX() const { return mX; }
  float getY() const { return mY; }
  float getZ() const { return mZ; }
  int getDetector() const { return mDetector; }
  int getPadRow() const { return mPadrow; }
  float getPadCol() const { return mPadcol; }
  float getTimeBin() const { return mTimebin; }

 protected:
  float mX, mY, mZ;
  int mDetector;
  int mPadrow;
  float mPadcol, mTimebin;

  static constexpr float xscale = 1.0 / (o2::trd::Geometry::cheight() + o2::trd::Geometry::cspace());
  static constexpr float xoffset = o2::trd::Geometry::getTime0(0);
  static constexpr float alphascale = 1.0 / o2::trd::Geometry::getAlpha();
};

ChamberSpacePoint::ChamberSpacePoint(o2::track::TrackParCov& trackpar, o2::trd::Geometry* geo)
  : mX(trackpar.getX()), mY(trackpar.getY()), mZ(trackpar.getZ())
{
  if (geo == nullptr)
    geo = o2::trd::Geometry::instance();

  int layer = int((mX - xoffset) * xscale);

  int sector = int(trackpar.getAlpha() * alphascale);
  while (sector < 0)
    sector += 18;
  while (sector >= 18)
    sector -= 18;

  int stack = geo->getStack(trackpar.getZ(), layer);
  if (stack < 0) {
    // if (draw)
    // printf("WARN: cannot determine stack for z = %f, layer = %i\n", trackpar.getZ(), layer);
    mDetector = -999;
    return;
  }

  mDetector = 30 * sector + 6 * stack + layer;
  auto pp = geo->getPadPlane(layer, stack);
  mPadrow = pp->getPadRowNumber(trackpar.getZ());
  mPadcol = pp->getPadColNumber(trackpar.getY());

  int rowMax = (stack == 2) ? 12 : 16;
  if (mPadrow < 0 || mPadrow >= rowMax) {
    // if (draw)
    printf("WARN: row  = %i for z = %f\n", mPadrow, trackpar.getZ());
    mDetector = -abs(mDetector);
  }

  // Mark points that are too far outside of the detector as questionable.
  if (mPadcol < -2.0 || mPadcol >= 145.0 ) {
    // if (draw)
    printf("WARN: col  = %f for y = %f\n", mPadcol, trackpar.getY());
    mDetector = -abs(mDetector);
  }

  // cout << "Created new ChamberSpacePoint with det=" << getDetector()
      //  << " row=" << getPadRow() << " col=" << getPadCol() << endl;
}

ostream& operator<<(ostream& os, const ChamberSpacePoint& p)
{
  os << "( " << setprecision(5) << p.getX() 
     << " / " << setprecision(5) << p.getY() 
     << " / " << setprecision(6) << p.getZ() << ") <-> ["
     << p.getDetector() << "." << p.getPadRow() 
     << " pad " << setprecision(5) << p.getPadCol() << "]";
  return os;
}



// ========================================================================
// ========================================================================
//
// RawDataSpan: container for TRD digits, tracklets, space points and
// other pieces of information to understand the TRD
//
// ========================================================================
// ========================================================================

// template <typename value_t, template <typename T> typename container_t>
template <typename value_t, typename container_t>
struct myrange {
  // typedef value_t value_type;
  // typedef typename container_t<value_t>::iterator iterator;
  typedef typename container_t::iterator iterator;

  iterator &begin() { return b; }
  iterator &end() { return e; }
  iterator b, e;

  size_t length() { return e - b; }
  void sort(bool (*comp)(const value_t &a, const value_t &b))
  {
    std::stable_sort(b, e, comp);
  }
};

class RawDataSpan
{
public:
  myrange<o2::trd::Digit, TTreeReaderArray<o2::trd::Digit>> digits;
  myrange<o2::trd::Tracklet64, TTreeReaderArray<o2::trd::Tracklet64>> tracklets;
  myrange<ChamberSpacePoint, std::vector<ChamberSpacePoint>> trackpoints;

  pair<int, int> getMaxADCsumAndChannel();
  int getMaxADCsum(){ return getMaxADCsumAndChannel().first; }

  int getDetector() { if (mDetector == -1) calculateCoordinates(); return mDetector; } 
  int getPadRow() { if (mPadRow == -1) calculateCoordinates(); return mPadRow; } 
  
protected:
  // The following variables cache the calculations of raw coordinates:
  //   non-negative values indicate the actual position
  //   -1 indicates that the values has not been calculated yet
  //   -2 indicates that the value is not unique in this span 
  int mDetector{-1};
  int mPadRow{-1};

  void calculateCoordinates();

};


void RawDataSpan::calculateCoordinates()
{
  for (auto& x : digits) {
    if (mDetector == -1) {
      mDetector = x.getDetector();
    } else if (mDetector != x.getDetector()) {
      mDetector = -2;
    }
    if (mPadRow == -1) {
      mPadRow = x.getPadRow();
    } else if (mPadRow != x.getPadRow()) {
      mPadRow = -2;
    }
  }
  for (auto& x : tracklets) {
    if (mDetector == -1) {
      mDetector = x.getDetector();
    } else if (mDetector != x.getDetector()) {
      mDetector = -2;
    }
    if (mPadRow == -1) {
      mPadRow = x.getPadRow();
    } else if (mPadRow != x.getPadRow()) {
      mPadRow = -2;
    }
  }
}



struct RawEvent : public RawDataSpan
{
  vector<o2::tpc::TrackTPC> tpctracks;
  vector<o2::dataformats::TrackTPCITS> tracks;
  vector<ChamberSpacePoint> evtrackpoints;
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
  os << sp.digits.length() << " digits, " 
     << sp.tracklets.length() << " tracklets, "
     << sp.trackpoints.length() << " track points";
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

bool order_spacepoint(const ChamberSpacePoint &a, const ChamberSpacePoint &b)
{
  if (a.getDetector() != b.getDetector())
    return a.getDetector() < b.getDetector();

  if (a.getPadRow() != b.getPadRow())
    return a.getPadRow() < b.getPadRow();

  if (a.getPadCol() != b.getPadCol())
    return a.getPadCol() < b.getPadCol();

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
  template<typename T> static uint32_t key(const T &x) { return x.getDetector(); }

  static bool comp_digits(const o2::trd::Digit &a, const o2::trd::Digit &b)
  { return key(a) != key(b); }

  static bool comp_tracklets(const o2::trd::Tracklet64 &a, const o2::trd::Tracklet64 &b)
  { return key(a) != key(b); }

  static bool comp_trackpoints(const ChamberSpacePoint &a, const ChamberSpacePoint &b)
  { return key(a) != key(b); }
};

struct ClassifierByPadRow
{
  template<typename T>
  static uint32_t key(const T &x) { return 100*x.getDetector() + x.getPadRow(); }

  static bool comp_digits(const o2::trd::Digit &a, const o2::trd::Digit &b)
  { return key(a) != key(b); }

  static bool comp_tracklets(const o2::trd::Tracklet64 &a, const o2::trd::Tracklet64 &b)
  { return key(a) != key(b); }

  static bool comp_trackpoints(const ChamberSpacePoint& a, const ChamberSpacePoint& b)
  { return key(a) != key(b); }
};

// struct ClassifierByMCM
// {
//   template<typename T>
//   static uint32_t key(const T &x) { 
//     return 1000*x.getDetector() + 10*x.getPadRow() + 4*(x.getROB()%2) + x.getMCM()%4;
//   }

//   template<>
//   static uint32_t key(const T &x) { 
//     return 1000*x.getDetector() + 10*x.getPadRow() + 4*(x.getROB()%2) + x.getMCM()%4;
//   }

//   static bool comp_digits(const o2::trd::Digit &a, const o2::trd::Digit &b)
//   { return key(a) != key(b); }

//   static bool comp_tracklets(const o2::trd::Tracklet64 &a, const o2::trd::Tracklet64 &b)
//   { return key(a) != key(b); }

//   static bool comp_trackpoints(const ChamberSpacePoint& a, const ChamberSpacePoint& b)
//   { return key(a) != key(b); }
// };

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
    event.trackpoints.sort(order_spacepoint);

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

    for (auto cur = event.trackpoints.b; cur != event.trackpoints.e; /* noop */) {
      auto nxt = std::adjacent_find(cur, event.trackpoints.e, classifier::comp_trackpoints);
      if (nxt != event.trackpoints.e) { ++nxt; }
      (*this)[classifier::key(*cur)].trackpoints.b = cur; 
      (*this)[classifier::key(*cur)].trackpoints.e = nxt;
      cur = nxt;
    }
    cout << "Found " << this->size() << " spans" << endl;
  }
};

// ========================================================================
// ========================================================================
//
// Track extrapolation
//
// ========================================================================
// ========================================================================

class TrackExtrapolator
{
 public:
  TrackExtrapolator();

  ChamberSpacePoint extrapolate(o2::track::TrackParCov& t, int layer);

 private:
  bool adjustSector(o2::track::TrackParCov& t);
  int getSector(float alpha);

  o2::trd::Geometry* mGeo;
  o2::base::Propagator* mProp;
};

TrackExtrapolator::TrackExtrapolator()
  // : mGeo(o2::trd::Geometry::instance())
  // : mProp(o2::base::Propagator::Instance())
{
  cout << "TrackExtrapolator: load geometry" << endl;
  o2::base::GeometryManager::loadGeometry();
  cout << "TrackExtrapolator: load B-field" << endl;
  o2::base::Propagator::initFieldFromGRP();

  cout << "TrackExtrapolator: instantiate geometry" << endl;

  mGeo = o2::trd::Geometry::instance();
  mGeo->createPadPlaneArray();
  mGeo->createClusterMatrixArray();

  mProp = o2::base::Propagator::Instance();
}

bool TrackExtrapolator::adjustSector(o2::track::TrackParCov& t)
{
  float alpha = mGeo->getAlpha();
  float xTmp = t.getX();
  float y = t.getY();
  float yMax = t.getX() * TMath::Tan(0.5f * alpha);
  float alphaCurr = t.getAlpha();
  if (fabs(y) > 2 * yMax) {
    printf("Skipping track crossing two sector boundaries\n");
    return false;
  }
  int nTries = 0;
  while (fabs(y) > yMax) {
    if (nTries >= 2) {
      printf("Skipping track after too many tries of rotation\n");
      return false;
    }
    int sign = (y > 0) ? 1 : -1;
    float alphaNew = alphaCurr + alpha * sign;
    if (alphaNew > TMath::Pi()) {
      alphaNew -= 2 * TMath::Pi();
    } else if (alphaNew < -TMath::Pi()) {
      alphaNew += 2 * TMath::Pi();
    }
    if (!t.rotate(alphaNew)) {
      return false;
    }
    if (!mProp->PropagateToXBxByBz(t, xTmp)) {
      return false;
    }
    y = t.getY();
    ++nTries;
  }
  return true;
}

int TrackExtrapolator::getSector(float alpha)
{
  if (alpha < 0) {
    alpha += 2.f * TMath::Pi();
  } else if (alpha >= 2.f * TMath::Pi()) {
    alpha -= 2.f * TMath::Pi();
  }
  return (int)(alpha * NSECTOR / (2.f * TMath::Pi()));
}

ChamberSpacePoint TrackExtrapolator::extrapolate(o2::track::TrackParCov& track, int layer)
{

  // ChamberSpacePoint pos;

  // if (draw)
  //   printf("Found ITS-TPC track with time %f us\n", trkITSTPC.getTimeMUS().getTimeStamp());

  double minPtTrack = 0.5;

  // const auto& track = trkITSTPC.getParamOut();
  if (fabs(track.getEta()) > 0.84 || track.getPt() < minPtTrack) {
    // no chance to find tracklets for these tracks
    return ChamberSpacePoint(-999);
  }

  if (!mProp->PropagateToXBxByBz(track, mGeo->getTime0(layer))) {
    // if (draw)
    printf("Failed track propagation into layer %i\n", layer);
    return ChamberSpacePoint(-999);
  }
  if (!adjustSector(track)) {
    // if (draw)
    printf("Failed track rotation in layer %i\n", layer);
    // return false;
  }

  // if (draw)
  // printf("Track has alpha of %f in layer %i. X(%f), Y(%f), Z(%f). Eta(%f), Pt(%f)\n", track.getAlpha(), layer, track.getX(), track.getY(), track.getZ(), track.getEta(), track.getPt());

  return ChamberSpacePoint(track, mGeo);

  // pos.x = track.getX();
  // pos.y = track.getY();
  // pos.z = track.getZ();

  // }
  //   ++nPointsTrigger;
  //   trackMap.insert(std::make_pair(std::make_tuple(geo->getDetector(layer, stack, sector), row, col), iTrack));
  //   int rowGlb = stack < 3 ? row + stack * 16 : row + 44 + (stack - 3) * 16; // pad row within whole sector
  //   int colGlb = col + sector * 144;                                         // pad column number from 0 to NSectors * 144
  //   hTracks[layer]->SetBinContent(rowGlb + 1, colGlb + 1, 4);
  // }

  // return true;
}

// ========================================================================
// ========================================================================
//
// the DataManager class
//
// ========================================================================
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
  RawEvent GetEvent();
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

  TrackExtrapolator extra;
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

    // Fix of position/slope should no longer be necessary
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

RawEvent DataManager::GetEvent()
{
  RawEvent ev;
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

  if (tracks) {
    for (auto &track : *tracks) {
      //   // auto tracktime = track.getTimeMUS().getTimeStamp();
      auto dtime = track.getTimeMUS().getTimeStamp() - evtime;
      if (dtime > mMatchTimeMinTPC && dtime < mMatchTimeMaxTPC) {
        ev.tracks.push_back(track);

        for(int ly=0; ly<6; ++ly) {
          auto point = extra.extrapolate(track.getParamOut(), ly);
          if (point.isValid()) {
            ev.evtrackpoints.push_back(point);
          }
        }
      }
    }
  }

  ev.trackpoints.b = ev.evtrackpoints.begin();
  ev.trackpoints.e = ev.evtrackpoints.end();

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
// ========================================================================
//
// Drawing routines
//
// ========================================================================
// ========================================================================

TVirtualPad *DrawPadRow(RawDataSpan &padrow, TVirtualPad *pad = NULL, TH2F *adcmap = NULL)
{
  // auto x = *padrow.digits.begin();
  // string desc = fmt::format("{:m}", x);
  string name = fmt::format("det{:03d}_row{:d}",
                            padrow.getDetector(), padrow.getPadRow());
  string desc = name;

  cout << "Plotting padrow " << name << endl;
  if (pad == NULL) {
    pad = new TCanvas(desc.c_str(), desc.c_str(), 1200, 500);
    pad->cd();
  }

  if (adcmap == NULL) {
    adcmap = new TH2F(name.c_str(), (desc + ";pad;time bin").c_str(), 144, 0., 144., 30, 0., 30.);
  }

  for (auto digit : padrow.digits) {
    if (digit.isSharedDigit()) { continue; }

    auto adc = digit.getADC();
    for (int tb = 0; tb < 30; ++tb) {
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

  TMarker mrk;
  mrk.SetMarkerStyle(20);
  mrk.SetMarkerSize(8);
  mrk.SetMarkerColor(kGreen);
  for (auto point : padrow.trackpoints)
  {
    mrk.DrawMarker(point.getPadCol(), 0);
    cout << point.getDetector() << " / "
         << point.getPadRow() << " / "
         << point.getPadCol() << endl;
    // auto pos = PadPosition(tracklet);
    // auto ypos = UncalibratedPad(tracklet);
    // auto slope = Slope(tracklet);
    // trkl.DrawLine(pos, 0, pos - 30*slope, 30);
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
