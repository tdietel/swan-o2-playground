
// #include "../macros/DataManager.C++"
#include "/Users/tom/alice/o2-playground/macros/DataManager.C++"

#include <string>
using namespace std;

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
    // if (padrow.digits.length() == 0 ) { continue; }
    // if (padrow.tracklets.length() == 0 ) { continue; }

    // int row = padrow.digits.begin()->getPadRow();
    int row = padrow.getPadRow();
    if (row < 0) {
      cout << "ERROR: span (" << padrow << ")" << " has pad row " << row << endl;
      continue;
    }
    pad[row]->cd();
    DrawPadRow(padrow, pad[row], adcdisp[row]);

    // txt.DrawTextNDC(0.2, 0.3, adcdisp[row]->GetName());
    txt.DrawTextNDC(0.15, 0.7, Form("row=%d", row));
    txt.DrawTextNDC(0.3, 0.7, Form("%zu tracklets", padrow.tracklets.length()));

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

// void draw_chamber(std::string dirname = "data/") 
void draw_chamber(std::string dirname = "/Users/tom/alice/data/alex-TRD-126-rec/")
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
      if (ev.digits.length() < 800000) { continue; }
      // if (ev.trackpoints.length() == 0) { continue; }

      cout << "============================================================================" << endl;
      cout << ev << endl;
      cout << "============================================================================" << endl;
      for (auto& point : ev.trackpoints) {
        cout << point << endl;
      }

      for (auto &[key, det] : RawDataPartitioner<ClassifierByDetector>(ev)) {
        if (det.digits.length() < 1) { continue; }
        // if (det.tracklets.length() < 1) { continue; }
        // if (det.trackpoints.length() < 1) {continue; }
        if (det.trackpoints.length() == 0) { continue; }

        // int detno = det.digits.begin()->getDetector();
        int detno = key;
        auto name = fmt::format(
          "det{:03d}_tf{}_ev{:04d}",
          detno, dman.GetTimeFrameNumber(), dman.GetEventNumber()
        );

        rawdisp->SetName(name.c_str());

        cout << endl << rawdisp->GetName() << endl;
        for (auto& point : det.trackpoints) {
          cout << point << endl;
        }

        for (auto& digit : det.digits) {
          cout << digit << endl;
        }

        DrawROC(det, rawdisp);
        // rawdisp->Modified();
        // rawdisp->Update();
        rawdisp->SaveAs((name+".pdf").c_str());

        cout << "--------------------------------------------------------------------------" << endl;

        // if (detno > 30) { return; }

        // sleep(10);

      }

      cout << "============================================================================" << endl;


      // return;


    } // event/trigger record loop
  } // time frame loop

  rawdispfile->Close();

}

