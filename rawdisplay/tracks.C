#include "../macros/DataManager.C++"

void tracks(std::string dirname = "/Users/tom/alice/data/alex-TRD-126-rec/")
{
  // ----------------------------------------------------------------------
  // instantiate the class that handles all the data access
  auto dman = DataManager(dirname);
  dman.SetMatchWindowTPC(-20., 30.);
  // auto dman = DataManager("/Users/tom/alice/data/ctf_evdisp/apass4/");

  // ----------------------------------------------------------------------
  // create output objects
  // TCanvas *rawdisp = new TCanvas();
  // rawdisp->Draw();

  TH2F* htime = new TH2F("htime", "htime;evtime;tracktime", 100, 0., 12000, 100, 0., 12000.);
  TH1F* hdtime = new TH1F("hdtime", "htime;tracktime - evtime;", 200, -150., 100.);
  TH1F* hdmtime = new TH1F("hdmtime", "htime (matched);tracktime - evtime;", 200, -150., 50.);
  TH1F* hditime = new TH1F("hditime", "hditime;tracktime - evtime;", 200, -150., 100.);

  // TNtuple *roi = new TNtuple("roi", "roi", "det:rob:mcm:ch:tbmax:nbhi:nblo:nbnb");
  // TNtuple *diginfo = new TNtuple("diginfo", "diginfo",
  //                                "det:rob:mcm:ch:adcsum:mean:rms:lmax");

  // ----------------------------------------------------------------------
  // loop over time frames and events
  while (dman.NextTimeFrame()) {
    int tfno = dman.GetTimeFrameNumber();

    while (dman.NextEvent()) {
      int evno = dman.GetEventNumber();
      auto ev = dman.GetEvent();

      auto tfid = dman.GetTimeFrameInfo();
      auto evtime = dman.GetTriggerTime();
      cout << dman.GetTriggerTime() << endl;

      // if (ev.digits.length() == 0) { continue; }
      // if (ev.digits.length() > 800000) { continue; }

      for (auto& track : *dman.GetTimeFrameTPCTracks()) {
        auto tracktime = track.getTime0()/5.0;
        htime->Fill(evtime, tracktime);
        hdtime->Fill(tracktime - evtime);
      }

      for (auto& track : *dman.GetTimeFrameTracks()) {
        auto tracktime = track.getTimeMUS().getTimeStamp();
        hditime->Fill(tracktime - evtime);
      }

      for (auto& track : ev.tpctracks) {
      //   // auto tracktime = track.getTimeMUS().getTimeStamp();
        auto tracktime = track.getTime0()/5.0;
        hdmtime->Fill(tracktime - evtime);
        //   //   cout << track << endl;
      }

      // for (auto &[key, mcm] : RawDataPartitioner<ClassifierByMCM>(ev)) {

      //   // Only consider MCMs without tracklets, but at least 3 digits.
      //   // if (mcm.digits.length() < 2) { continue; }
      //   // if (mcm.tracklets.length() >= 1) { continue; }


      // }
      // roi->Write();
      // diginfo->Write();

    } // event/trigger record loop
  }   // time frame loop

  hditime->SetLineColor(kRed);

  hditime->DrawClone();
  hdtime->DrawClone("same");

  TFile *outfile = new TFile("trd_track_matches.root", "RECREATE");
  htime->Write();
  hdtime->Write();
  hditime->Write();
  hdmtime->Write();
  outfile->Close();
}
