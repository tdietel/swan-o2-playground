
#include "../macros/DataManager.C+"

TVirtualPad *DrawROC(RawDataSpan &det, TPad *canvas = NULL)
{
  int height = 100;
  int width = 800;
  float bottom = 0.5;
  float top = 0.2;

  auto x = *det.digits.begin();
  string desc = fmt::format("Detector {:03d}", x.getDetector());
  string name = fmt::format("det{:03d}", x.getDetector());

  if (canvas == NULL) {
    canvas = new TCanvas(desc.c_str(), desc.c_str(), width, 16*height+top+bottom);
  } else {
    canvas->SetName(name.c_str());
    canvas->SetTitle(desc.c_str());
  }

  TPad* pad[16];
  TH2F* adcdisp[16];
  
  for (int i=0; i<16; ++i) {
    // Calculate vertical pad position
    float hunit = 1.0 / (16.0 + top + bottom);
    float lo = i == 0 ? 0.0 : (float(i) + bottom) * hunit;
    float hi = i == 15 ? 1.0 : (float(i+1) + bottom) * hunit;

    name = fmt::format("padrow_{:02}", i);
    canvas->cd(0);
    pad[i] = (TPad *)gROOT->FindObject(name.c_str());
    if (pad[i] == NULL) {
      cout << "Create " << name << endl;
      pad[i] = new TPad(name.c_str(), "", 0, lo, 1.0, hi);
      pad[i]->SetTopMargin(i == 15 ? top / (1.0+top) : 0.02);
      pad[i]->SetBottomMargin(i==0 ? bottom / (1.0+bottom) : 0.02);
      pad[i]->SetRightMargin(0.04);
      pad[i]->SetLogz();
    }
    pad[i]->Draw();
    pad[i]->cd();

    name = fmt::format("adcdisp_padrow{:02}", i);
    adcdisp[i] = (TH2F *)gROOT->FindObject(name.c_str());
    if (adcdisp[i] == NULL) {
      cout << "Create " << name << endl;
      adcdisp[i] = new TH2F(name.c_str(), ";;time bin", 144, 0., 144., 30, 0., 30.);
      adcdisp[i]->SetStats(0);

      adcdisp[i]->GetXaxis()->SetLabelFont(43);
      adcdisp[i]->GetXaxis()->SetLabelSize(14);
      adcdisp[i]->GetXaxis()->SetTitleFont(43);
      adcdisp[i]->GetXaxis()->SetTitleSize(16);
      adcdisp[i]->GetXaxis()->SetTitleOffset(1.4);

      adcdisp[i]->GetYaxis()->SetTickSize(0.00);
      adcdisp[i]->GetYaxis()->SetTicks("-");
      adcdisp[i]->GetYaxis()->SetLabelFont(43);
      adcdisp[i]->GetYaxis()->SetLabelSize(16);
    } 
      
    adcdisp[i]->Reset();
    adcdisp[i]->SetMinimum(0.9);
    adcdisp[i]->SetMaximum(1024.);
    adcdisp[i]->Draw("colz");
  }

  // for (int i = 0; i < 16; ++i)
  // {
  // }

  TText txt;
  txt.SetTextFont(43);
  txt.SetTextSize(18);
  for (auto &[key, padrow] : RawDataPartitioner<ClassifierByPadRow>(det))
  {
    if (padrow.digits.length() == 0 ) { continue; }
    // if (padrow.tracklets.length() == 0 ) { continue; }

    int row = padrow.digits.begin()->getPadRow();
    pad[row]->cd();
    DrawPadRow(padrow, pad[row], adcdisp[row]);

    // txt.DrawTextNDC(0.2, 0.3, adcdisp[row]->GetName());
    txt.DrawTextNDC(0.15, 0.7, Form("row=%d", row));
    txt.DrawTextNDC(0.3, 0.7, Form("%d tracklets", padrow.tracklets.length()));

    // txt.DrawTextNDC(0.5, 0.3, pad[row]->GetName());
    // txt.DrawTextNDC(0.2, 0.6, Form("row=%d", row));
    // auto c = DrawPadRow(padrow);
    // c->Update();
    // c->Draw();

  }

  // for (int i=1; i<=16; ++i) {
  //   cout << i << " -> " << canvas->GetPad(i) << endl;
  //   // auto iter = pad->GetPad(i)->GetListOfPrimitives()->MakeIterator();
  //   // while(auto prim = iter->Next()) {
  //   //   cout << prim->GetName() << endl;
  //   // }
  // }

  return canvas;
}

void draw_chamber(std::string dirname="data/")
{

  // ----------------------------------------------------------------------
  // instantiate the class that handles all the data access
  // auto dman = DataManager("./");
  auto dman = DataManager(dirname);

  // ----------------------------------------------------------------------
  // create output objects

  TCanvas* rawdisp = new TCanvas("bar","blub", 900, 1500);
  TFile* rawdispfile = new TFile("rawdisplay.root", "RECREATE");


  // rawdisp->Update();
  // rawdisp->Draw();
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

      for (auto &[key, det] : RawDataPartitioner<ClassifierByDetector>(ev)) {
        if (det.digits.length() < 3) { continue; }
        // if (det.tracklets.length() >= 1) { continue; }

        int detno = det.digits.begin()->getDetector();
        rawdisp->SetName(fmt::format(
          "det{:03d}_tf{}_ev{:04d}",
          detno, dman.GetTimeFrameNumber(), dman.GetEventNumber()
        ).c_str());
                         
        DrawROC(det, rawdisp);
        rawdisp->SaveAs(".pdf");

        // if (detno > 30) { return; }

      }

      return;

      // rawdisp->Update();
      // rawdisp->Draw();
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

}

