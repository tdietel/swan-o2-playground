#include "/Users/tom/alice/o2-playground/macros/DataManager.C++"

#include <iostream>
#include <TH1.h>
#include <TH2.h>

using namespace std;

void tracks(std::string dirname = "/Users/tom/alice/data/alex-TRD-126-rec/")
{
  // ----------------------------------------------------------------------
  // instantiate the class that handles all the data access
  cout << "Creating data manager" << endl;
  auto dman = DataManager(dirname);
  dman.SetMatchWindowTPC(-20., 30.);
  // auto dman = DataManager("/Users/tom/alice/data/ctf_evdisp/apass4/");

  // ----------------------------------------------------------------------
  // create output objects
  // TCanvas *rawdisp = new TCanvas();
  // rawdisp->Draw();

  TFile *outfile = new TFile("trd_track_matches.root", "RECREATE");

  TH2F* htime = new TH2F("htime", "htime;evtime;tracktime", 100, 0., 12000, 100, 0., 12000.);
  TH1F* hdtime = new TH1F("hdtime", "htime;tracktime - evtime;", 200, -150., 100.);
  TH1F* hdmtime = new TH1F("hdmtime", "htime (matched);tracktime - evtime;", 200, -150., 50.);
  TH1F* hditime = new TH1F("hditime", "hditime;tracktime - evtime;", 200, -150., 100.);

  TH1F* hpadres = new TH1F("hpadres", "pad resolution;col_{track}-col_{tracklet};N", 100, -20., 20.);
  TH2F* hpadressm = new TH2F("hpadressm", "pad resolution by SM;col_{track}-col_{tracklet};sector", 100, -20., 20., 18, -0.5, 17.5);
  TH2F* hpadpos = new TH2F("hpadpos", "pad resolution;col_{track};col_{tracklet}", 288, 0., 144., 288, 0., 144.);

  TNtuple* res = new TNtuple("res","res", "det:sm:trow:tcol:prow:pcol");

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

      cout << ev << endl;

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

      // for (auto& point : ev.trackpoints) {
      //   cout << point << endl;
      // }
      for (auto &[key, row] : RawDataPartitioner<ClassifierByDetector>(ev)) {

        for (auto& point : row.trackpoints) {
          for (auto& tracklet : row.tracklets) {
            float track_pos = point.getPadCol();
            float tracklet_pos = PadPosition(tracklet);

            res->Fill(tracklet.getHCID()/2, tracklet.getHCID()/60, 
              tracklet.getPadRow(), PadPosition(tracklet), point.getPadRow(), point.getPadCol());
            hpadres->Fill(track_pos - tracklet_pos);
            hpadressm->Fill(track_pos - tracklet_pos, tracklet.getHCID()/60);
            hpadpos->Fill(track_pos, tracklet_pos);
            if ( fabs(track_pos - tracklet_pos) < 10) {
              cout << point << endl << tracklet << endl;
            }
          }
        }
      }

      //   cout << point << endl;
      // }


        //   // Only consider MCMs without tracklets, but at least 3 digits.
        //   // if (mcm.digits.length() < 2) { continue; }
        //   // if (mcm.tracklets.length() >= 1) { continue; }

        // }
        // roi->Write();
        // diginfo->Write();

    } // event/trigger record loop
  }   // time frame loop

  hditime->SetLineColor(kRed);

  // hditime->DrawClone();
  // hdtime->DrawClone("same");

  hpadressm->Draw("colz");

  // TFile *outfile = new TFile("trd_track_matches.root", "RECREATE");
  htime->Write();
  hdtime->Write();
  hditime->Write();
  hdmtime->Write();
  hpadres->Write();
  hpadressm->Write();
  hpadpos->Write();
  res->Write();
  outfile->Close();
}
