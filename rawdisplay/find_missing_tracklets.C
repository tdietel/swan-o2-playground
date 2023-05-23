#include "../macros/DataManager.C+"

void find_missing_tracklets(std::string dirname="data/")
{
  // ----------------------------------------------------------------------
  // instantiate the class that handles all the data access
  auto dman = DataManager(dirname);
  // auto dman = DataManager("/Users/tom/alice/data/ctf_evdisp/apass4/");

  // ----------------------------------------------------------------------
  // create output objects
  TCanvas* rawdisp = new TCanvas();
  TFile* rawdispfile = new TFile("rawdisplay.root", "RECREATE");

  // ----------------------------------------------------------------------
  // loop over time frames and events
  while ( dman.NextTimeFrame() ) {
    int tfno = dman.GetTimeFrameNumber();

    while (dman.NextEvent()) {
      int evno = dman.GetEventNumber();
      auto ev = dman.GetEvent();
      if (ev.digits.length() == 0) { continue; }
      if (ev.digits.length() > 800000) { continue; }

      for (auto &[key, mcm] : RawDataPartitioner<ClassifierByMCM>(ev)) {

        // Only consider MCMs without tracklets, but at least 3 digits.
        if (mcm.digits.length() < 3) { continue; }
        if (mcm.tracklets.length() >= 1) { continue; }

        // ignore MCMs where the ADC sum of the highest digit is too small
        if (mcm.getMaxADCsum() < 350) continue;

        auto pad = DrawMCM(mcm, rawdisp);
        // pad->SaveAs(".pdf");
        auto name = string(pad->GetTitle());
        name += fmt::format("  tf {}  ev {}", tfno, evno);

        rawdispfile->WriteObject(pad, name.c_str() );
      }

    } // event/trigger record loop
  } // time frame loop

  rawdispfile->Close();

}

