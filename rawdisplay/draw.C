
#include "../macros/DataManager.C++"

// void draw(std::string dirname="data/")
void draw(std::string dirname = "/Users/tom/alice/data/alex-TRD-126-rec/")
{

  // ----------------------------------------------------------------------
  // instantiate the class that handles all the data access
  // auto dman = DataManager("./");
  auto dman = DataManager(dirname);
  // auto dman = DataManager("/Users/tom/alice/data/ctf_evdisp/apass4/");

  // ----------------------------------------------------------------------
  // create output objects
  TH2F* poscalib = new TH2F("poscalib", "poscalib", 21, -0.5, 20.5, 60, -5., 25.);
  TCanvas* cnv = new TCanvas("cnv_padrow", "cnv_padrow", 800,600);

  TCanvas* rawdisp = new TCanvas();
  TFile* rawdispfile = new TFile("rawdisplay.root", "RECREATE");

  // SpacePointConverter conv;

  // ----------------------------------------------------------------------
  // loop over time frames and events
  // comparator_padrow comparator;
  // comparator_det comparator;

  int ndrawn = 0;

  while ( dman.NextTimeFrame() ) {
    while ( dman.NextEvent() ) {

      auto ev = dman.GetEvent();
      if (ev.digits.length() == 0) { continue; }
      if (ev.digits.length() > 800000) { continue; }

      // for (auto det : RawDataPartitionByDetector(dman.GetEvent()))
      // for (auto &[key, det] : RawDataPartitionByDetector(dman.GetEvent()))
      // for (auto &[key, det] : RawDataPartitioner<ClassifierByPadRow>(dman.GetEvent()))
      for (auto &[key, det] : RawDataPartitioner<ClassifierByMCM>(dman.GetEvent()))
      {
        if (det.digits.length() < 3) { continue; }
        if (det.tracklets.length() < 1) { continue; }

        bool accept_roc = false;
        for (auto tracklet : det.tracklets)
        {
          if (PadPositionMCM(tracklet) > 3 && PadPositionMCM(tracklet) < 18)
          {
            accept_roc = true;
          }
        }
        // if (accept_roc == false) { continue; }

        cout << "=====================" << endl;
        cout << det << endl;
        int maxch = -1, maxval = -1;
        for (auto digit : det.digits)
        {
          // if (digit.getADCsum() < 400) { continue; }
          // if (digit.getChannel() <= 1 || digit.getChannel() >= 19) { continue; }
          if (digit.getADCsum() > maxval)
          {
            maxval = digit.getADCsum();
            maxch = digit.getChannel();
          }
          cout << digit << endl;
        }

        for (auto tracklet : det.tracklets)
        {
          // cout << tracklet << endl;
          poscalib->Fill(maxch, PadPosition(tracklet));
          fmt::print("{:c}\n", tracklet);
          cout << "Tracklet Det " << tracklet.getDetector()
               << " ROB " << tracklet.getROB()
               << " MCM " << tracklet.getMCM()
               << " pos=" << tracklet.getPosition() 
               << "=" << PadPositionMCM(tracklet)
               << " slope=" << tracklet.getSlope()
               << " y=" << tracklet.getUncalibratedY()
               << " dy=" << tracklet.getUncalibratedDy(15.0)
               << endl;
        }

        auto pad = DrawMCM(det, rawdisp);
        // pad->SaveAs(".pdf");
        auto name = string(pad->GetTitle());
        name += fmt::format("  tf {}  ev {}", 
                            dman.GetTimeFrameNumber(), dman.GetEventNumber());

        rawdispfile->WriteObject(pad, name.c_str() );
        // if (ndrawn < 50) {
        //   pad->Draw();
        //   ++ndrawn;
        // }
      }

      // while ( dman.NextMCM() ) {

      //   cout << "-----" << endl;
      //   for (auto& digit: dman.Digits()) {
      //     // if (digit.getADCsum() < 400) { continue; }
      //     // if (digit.getChannel() <= 1 || digit.getChannel()>=19) {continue;}

      //     cout << digit << endl;
      //   }
      // }
      // for (auto& seq : dman.DigitsByPadRow()) {

        // padrow->Reset();
        // padrow->SetStats(0);
        // padrow->GetXaxis()->SetRange(seq.begin()->getPadCol(),
        // (seq.end()-1)->getPadCol()+2);


        // int det = seq.begin()->getDetector();
        // int row = seq.begin()->getPadRow();

        // cout << (*seq.begin()) << " digits: "
        //      << (seq.end()-seq.begin()) << endl;

        // for (auto& dig : seq) {
        //   cout << "   " << dig << endl;

        //   auto adc = dig.getADC();
        //   for (int i=0;i<30;i++) {
        //     padrow->Fill(dig.getPadCol(), i, adc[i]);
        //   }
        // }

        // // new TCanvas();
        // // auto p = padrow->Clone(Form("padrow_%03d_%02d",det,row));
        // //p->Draw("colz");
        // // p->Write();


        // padrow->SetTitle(Form("Det %03d row %02d;pad;time bin",det,row));
        // padrow->Draw("colz");
        // padrow->Draw("text,same");


        // TMarker m;
        // m.SetMarkerStyle(20);

        // for (auto& hit : dman.Hits()) {

        //   // only use hits in current detector
        //   if (hit.GetDetectorID()!=det) continue;

        //   // convert xyz to pad row/col/timebin coordinates
        //   auto rct = conv.Hit2RowColTime(hit);

        //   // restrict to current padrow
        //   if (rct[0]<float(row) || rct[0]>float(row+1)) continue;

        //   cout << hit << "   " << rct[0] << ":" << rct[1] << ":" << rct[2] << endl;
        //   m.DrawMarker(rct[1], rct[2]);
        // }

        // cnv->SaveAs(Form("padrow_%03d_%02d.pdf",det,row));

      // } // padrow loop


      // return;

    } // event/trigger record loop
  } // time frame loop

  rawdispfile->Close();

  cnv->cd();
  poscalib->Draw("colz");
  TLine l;
  l.DrawLine(0,0,20,20);

}


// TPad* DrawMCM(RawDataSpan& mcm, TPad* pad)
// {
//   auto x = *mcm.digits.begin();
//   string desc = fmt::format("{:m}", x);
//   string name = fmt::format("det{:03d}_rob{:d}_mcm{:02d}",
//                             x.getDetector(), x.getROB(), x.getMCM());

//   if (pad == NULL)
//   {
//     pad = new TCanvas(desc.c_str(), desc.c_str(), 800, 600);
//   } else {
//     pad->SetName(name.c_str());
//     pad->SetTitle(desc.c_str());
//   }
//   pad->cd();

//   TH2F *digit_disp = new TH2F(desc.c_str(), (desc + ";ADC channel;time bin").c_str(), 21, 0., 21., 30, 0., 30.);

//   for (auto digit : mcm.digits) {
//     auto adc = digit.getADC();
//     for (int tb=0; tb<30; ++tb) {
//       digit_disp->Fill(digit.getChannel(), tb, adc[tb]);
//     }
//   }
//   digit_disp->SetStats(0);
//   digit_disp->Draw("colz");

//   TLine trkl;
//   trkl.SetLineColor(kRed);
//   trkl.SetLineWidth(3);

//   for (auto tracklet : mcm.tracklets)
//   {
//     auto pos = PadPositionMCM(tracklet);
//     auto slope = Slope(tracklet);
//     trkl.DrawLine(pos,0, pos + 30*slope, 30);
//   }
  
//   return pad;
// }

// void DrawROC(RawDataSpan roc)
// {}
