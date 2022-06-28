
#include "../macros/DataManager.C"
#include "../macros/HelperFunctions.C"

#include <TPad.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TLine.h>


TPad* DrawPadRow(RawDataSpan &padrow, TPad* pad=NULL)
{
  auto x = *padrow.digits.begin();
  string desc = fmt::format("{:m}", x);
  string name = fmt::format("det{:03d}_rob{:d}_mcm{:02d}",
                            x.getDetector(), x.getROB(), x.getMCM());

  if (pad == NULL) {
    pad = new TCanvas(desc.c_str(), desc.c_str(), 1200, 500);
  }
  else {
    pad->SetName(name.c_str());
    pad->SetTitle(desc.c_str());
  }
  pad->cd();

  TH2F *adcmap = new TH2F(name.c_str(), (desc + ";pad;time bin").c_str(), 144, 0., 144., 30, 0., 30.);

  for (auto digit : padrow.digits) {
    if (digit.isSharedDigit()) {
      continue;
    }

    auto adc = digit.getADC();
    for (int tb = 0; tb < 30; ++tb) {
      adcmap->Fill(digit.getPadCol(), tb, adc[tb]);
    }
  }
  adcmap->SetStats(0);
  adcmap->Draw("colz");

  TLine trkl;
  trkl.SetLineColor(kGreen + 2);

  for (auto tracklet : padrow.tracklets) {
    auto pos = PadPosition(tracklet);
    auto slope = Slope(tracklet);
    trkl.DrawLine(pos, 0, pos + slope, 30);
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
