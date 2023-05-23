#include "../macros/DataManager.C+"

pair<float,float> AdcMeanRms(const o2::trd::Digit& digit)
{
  auto adc = digit.getADC();

  int sum = 0, sumsq = 0, n = 0;
  for(int tb=0; tb<30;tb++) {
    sum += adc[tb];
    sumsq += adc[tb] * adc[tb];
    if (adc[tb]) ++n;
  }

  float mean = float(sum) / float(n);
  float rms = float(sumsq) / float(n) - mean*mean;
  return make_pair(mean, rms);
}


void pulseheight(std::string dirname="./")
{
  // ----------------------------------------------------------------------
  // instantiate the class that handles all the data access
  auto dman = DataManager(dirname);
  // auto dman = DataManager("/Users/tom/alice/data/ctf_evdisp/apass4/");

  // ----------------------------------------------------------------------
  // create output objects
  TCanvas* rawdisp = new TCanvas();
  TFile* outfile = new TFile("pulseheight.root", "RECREATE");
  TNtuple* roi = new TNtuple("roi", "roi", "det:rob:mcm:ch:tbmax:nbhi:nblo:nbnb");
  TNtuple* diginfo = new TNtuple("diginfo", "diginfo",      
    "det:rob:mcm:ch:adcsum:mean:rms:lmax");

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
        if (mcm.digits.length() < 2) { continue; }
        // if (mcm.tracklets.length() >= 1) { continue; }

        int adcsum[21];
        bool localmax[21];
        for (int ch=0; ch<21; ++ch) {
          adcsum[ch] = 0;
          localmax[ch] = false;
        }

        for ( auto& digit : mcm.digits ) {
          int ch = digit.getChannel();
          adcsum[ch] = digit.getADCsum();
          if (ch > 0 && adcsum[ch] > adcsum[ch - 1]) {
            localmax[ch] = true;
            localmax[ch-1] = false;
          }
        }
        localmax[20] = false;

        int tbmax, nbhi, bhlo, nbnb;
        for ( auto& digit : mcm.digits ) {
          int ch = digit.getChannel();
          auto meanrms = AdcMeanRms(digit);
          diginfo->Fill(digit.getDetector(), digit.getROB(), digit.getMCM(),
                         digit.getChannel(),
                         adcsum[ch], meanrms.first, meanrms.second, localmax[ch]);
        }

        for ( int ch = 1; ch < 20; ++ch) {
          if ( ! localmax[ch] ) continue;
          if (adcsum[ch+1] > adcsum[ch-1]) {
            int nbnb = ch<19 ? adcsum[ch+2] : 0;
            roi->Fill(0,0,0,ch,adcsum[ch], adcsum[ch+1], adcsum[ch-1], nbnb);
          } else {
            int nbnb = ch>1 ? adcsum[ch-2] : 0;
            roi->Fill(0,0,0,ch, adcsum[ch], adcsum[ch-1], adcsum[ch+1], nbnb);
          }
        }

          // ignore MCMs where the ADC sum of the highest digit is too small
          // if (mcm.getMaxADCsum() < 350) continue;

          // auto pad = DrawMCM(mcm, rawdisp);
          // pad->SaveAs(".pdf");
          // auto name = string(pad->GetTitle());
          // name += fmt::format("  tf {}  ev {}", tfno, evno);

          // outfile->WriteObject(pad, name.c_str() );
      }
      roi->Write();
      diginfo->Write();

    } // event/trigger record loop
  } // time frame loop

  outfile->Close();

}

