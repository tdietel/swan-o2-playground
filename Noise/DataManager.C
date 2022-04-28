
// #include <TPython.h>

int DigitIdx(const o2::trd::Digit& a)
{
  int ret = a.getDetector();
  ret *= 1000; ret += a.getPadRow();
  ret *= 1000; ret += a.getPadCol();
  return ret;
}

bool DigitCompare(const o2::trd::Digit& a, const o2::trd::Digit& b)
{
  return DigitIdx(a) < DigitIdx(b);
}


class DataManager
{

public:

  DataManager(std::string dir="data/")
  : mainfile(0), datatree(0), datareader(0),
    hits(0), digits(0), tracklets(0), trgrecords(0),
    tfno(0), evno(0), iEvent(-1)
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
    
    // set up the branches we want to read
    digits     = new TTreeReaderArray<o2::trd::Digit>(*datareader, "TRDDigit");
    tracklets  = new TTreeReaderArray<o2::trd::Tracklet64>(*datareader, "Tracklet");
    trgrecords = new TTreeReaderArray<o2::trd::TriggerRecord>(*datareader, "TrackTrg");

    // ConnectMCHitsFile(dir+"o2sim_HitsTRD.root");

  }

  bool NextTimeFrame()
  {
    if (datareader->Next() ) {

      cout << "## Time frame " << tfno << endl;

      iEvent = 0;
      evno = 0;
      tfno++;
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
    //   cout << "no hits found for event " << tfno << ":" << evno << endl;
    //   return false;
    // }

    auto evrec = trgrecords->At(iEvent);
    cout << evrec << endl;

    cout << "## Event " << tfno << ":" << evno << ":  "
        //  << hits->GetSize() << " hits   "
         << evrec.getNumberOfDigits() << " digits" << endl;

    event_digits.b = digits->begin() + evrec.getFirstDigit();
    event_digits.e = event_digits.begin() + evrec.getNumberOfDigits();
    std::stable_sort(event_digits.b,event_digits.e,DigitCompare);

    // event_hits.b  = hits->begin();
    // event_hits.e  = hits->end();

    evno++;
    iEvent++;
    return true;
  }

  // o2::trd::TriggerRecord GetTriggerRecord() {return trgrecords->At(iEvent);}



  // A struct to be used for range-based for loops
  // A template would be nice, but I did not manage to make it work.
  struct DigitRange {
    TTreeReaderArray<o2::trd::Digit>::iterator& begin() { return b; }
    TTreeReaderArray<o2::trd::Digit>::iterator& end() { return e; }
    TTreeReaderArray<o2::trd::Digit>::iterator b, e, c;

    /// I had the idea to implement the iterator in C++, but didn't manage
    // DigitRange* __iter__() { c=b; return this; }
    // o2::trd::Digit __next__() {
    //   if (c==e) {
    //     TPython.Eval("raise StopIteration");
    //   } else {
    //     return *(c++);
    //   }
    // }
  };


  DigitRange event_digits, padrow_digits;
  DigitRange Digits() {return event_digits;}

  // A struct to be used for range-based for loops
  // A template would be nice, but I did not manage to make it work.
  struct HitRange {
    TTreeReaderArray<o2::trd::Hit>::iterator& begin() { return b; }
    TTreeReaderArray<o2::trd::Hit>::iterator& end() { return e; }
    TTreeReaderArray<o2::trd::Hit>::iterator b, e;
   };

   HitRange event_hits;
   HitRange Hits() {return event_hits;}


protected:
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

  // void ConnectDigitsFile(std::string fname)
  // {

  //   fdigits = new TFile(fname.c_str());
  //   fdigits->GetObject("o2sim", trdigits);

  //   rdrdigits = new TTreeReader(trdigits);

  //   // set up the branches we want to read
  //   digits = new TTreeReaderArray<o2::trd::Digit>(*rdrdigits, "TRDDigit");
  //   trgrecords = new TTreeReaderArray<o2::trd::TriggerRecord>(*rdrdigits, "TriggerRecord");
  // }

private:
  // std::list<TFile*> files;
  TFile* mainfile;
  // TFile *fhits, *fdigits, *ftracklets;
  // TTree *trhits, *trdigits, *trtracklets, *trtrgrec;
  TTree* datatree; // tree and friends from digits, tracklets files
  // TTreeReader *rdrhits, *rdrdigits, *rdrtracklets, *rdrtrgreg;
  TTreeReader *datareader;

  TTreeReaderArray<o2::trd::Hit>* hits;
  TTreeReaderArray<o2::trd::Digit>* digits;
  TTreeReaderArray<o2::trd::Tracklet64>* tracklets;
  TTreeReaderArray<o2::trd::TriggerRecord>* trgrecords;

  int tfno, evno, iEvent;
};
