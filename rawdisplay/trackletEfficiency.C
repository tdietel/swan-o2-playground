// ============================================================================
// Tracklet efficiency determination
//
// File written by Ole Schmidt, included here as reference for extrapolating
// tracks and matching them to digits.
// ============================================================================

// #if !defined(__CLING__) || defined(__ROOTCLING__)
// ROOT header
#include <TROOT.h>
#include <TChain.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLatex.h>
#include <TMath.h>
// O2 header
#include "DataFormatsTRD/TriggerRecord.h"
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/Constants.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTRD/NoiseCalibration.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "TRDBase/Geometry.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonDataFormat/TFIDInfo.h"
#include "CommonDataFormat/InteractionRecord.h"
// #include "TRDQC/StatusHelper.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonUtils/TreeStream.h"
#include "CommonUtils/TreeStreamRedirector.h"
// #include "CommonUtils/ConfigurableParam.h"
// #include "CommonUtils/ConfigurableParamHelper.h"

#include <set>
#include <cmath>
#include <unordered_map>
#include <algorithm>
// #endif

using namespace o2::trd;
using namespace o2::trd::constants;

// some settings
int searchRoadCol = 3;
int searchRoadRow = 1;
float minPtTrack = .5f;
int minAdcCount = 0;
bool doXOR = true;



TFile* fInTrklt = nullptr;
TFile* fInDigit = nullptr;
TFile* fInTracks = nullptr;
TFile* fInTFID = nullptr;
TTree* treeTrklt = nullptr;
TTree* treeDigit = nullptr;
TTree* treeTrks = nullptr;
o2::trd::Geometry* geo = nullptr;
std::vector<o2::dataformats::TrackTPCITS> tracks, *tracksInPtr{&tracks};
std::vector<TriggerRecord> trigIn, *trigInPtr{&trigIn};
std::vector<Tracklet64> tracklets, *trackletsInPtr{&tracklets};
std::vector<Digit> digitIn, *digitInPtr{&digitIn};
std::vector<o2::dataformats::TFIDInfo>* tfids = nullptr;

std::vector<TH2F*> hDigits;
std::vector<TH2F*> hTracklets;
std::vector<TH2F*> hTracks;

TCanvas* c = nullptr;
TLine* l = nullptr;

bool initDone = false;
bool initDrawingDone = false;
size_t lastEntry = 0;
size_t lastTrigger = -1;

bool prepareInput()
{
  fInTrklt = TFile::Open("trdtracklets.root");
  treeTrklt = (TTree*) fInTrklt->Get("o2sim");
  treeTrklt->SetBranchAddress("TrackTrg", &trigInPtr);
  treeTrklt->SetBranchAddress("Tracklet", &trackletsInPtr);

  fInDigit = TFile::Open("trddigits.root");
  treeDigit = (TTree*) fInDigit->Get("o2sim");
  treeDigit->SetBranchAddress("TRDDigit", &digitInPtr);

  fInTracks = TFile::Open("o2match_itstpc.root");
  treeTrks = (TTree*) fInTracks->Get("matchTPCITS");
  treeTrks->SetBranchAddress("TPCITS", &tracksInPtr);

  fInTFID = TFile::Open("o2_tfidinfo.root");
  if (fInTFID) {
    // for a simulation this file is not available
    tfids = (std::vector<o2::dataformats::TFIDInfo>*) fInTFID->Get("tfidinfo");
  }

  o2::base::GeometryManager::loadGeometry();
  o2::base::Propagator::initFieldFromGRP();

  geo = o2::trd::Geometry::instance();
  geo->createPadPlaneArray();
  geo->createClusterMatrixArray();

  if (treeDigit->GetEntries() != treeTrklt->GetEntries() || treeTrklt->GetEntries() != treeTrks->GetEntries()) {
    printf("Number of tree entries is not equal\n");
    return false;
  }
  if (tfids && tfids->size() != static_cast<size_t>(treeDigit->GetEntries())) {
    printf("We don't have TFIDs (%lu) for every TF (%lld) available\n", tfids->size(), treeDigit->GetEntries());
    return false;
  }
  return true;
}

void drawTrdGrid(TLine* line)
{
  line->DrawLine(15.5, 0, 15.5, 2592);
  line->DrawLine(31.5, 0, 31.5, 2592);
  line->DrawLine(43.5, 0, 43.5, 2592);
  line->DrawLine(59.5, 0, 59.5, 2592);
  for (int iSec = 1; iSec < 18; ++iSec) {
    float yPos = iSec * 144 - 0.5;
    line->DrawLine(0, yPos, 76, yPos);
  }
}

void createTrdPadHistsPerLayer()
{
  for (int iLayer = 0; iLayer < 6; ++iLayer) {
    hTracklets.push_back(new TH2F(Form("tracklets%i", iLayer), Form("Tracklet count per pad in layer %i;stack;sector", iLayer), 76, -0.5, 75.5, 2592, -0.5, 2591.5));
    hTracks.push_back(new TH2F(Form("tracks%i", iLayer), Form("Track count per pad in layer %i;stack;sector", iLayer), 76, -0.5, 75.5, 2592, -0.5, 2591.5));
    hDigits.push_back(new TH2F(Form("digits%i", iLayer), Form("Digit count per pad in layer %i;stack;sector", iLayer), 76, -0.5, 75.5, 2592, -0.5, 2591.5));
    auto xax = hDigits.back()->GetXaxis();
    xax->SetBinLabel(8, "0");
    xax->SetBinLabel(24, "1");
    xax->SetBinLabel(38, "2");
    xax->SetBinLabel(52, "3");
    xax->SetBinLabel(68, "4");
    xax->SetTicks("-");
    xax->SetTickSize(0.01);
    xax->SetLabelSize(0.045);
    xax->SetLabelOffset(0.01);
    xax->SetTitleOffset(-1);
    auto yax = hDigits.back()->GetYaxis();
    for (int iSec = 0; iSec < 18; ++iSec) {
      auto lbl = std::to_string(iSec);
      yax->SetBinLabel(iSec * 144 + 72, lbl.c_str());
    }
    yax->SetTicks("-");
    yax->SetTickSize(0.01);
    yax->SetLabelSize(0.045);
    yax->SetLabelOffset(0.01);
    yax->SetTitleOffset(1.4);
    hDigits.back()->SetStats(0);
    hTracklets.back()->SetStats(0);
    hTracks.back()->SetStats(0);
  }
}

void drawResults()
{
  if (!initDrawingDone) {
    c = new TCanvas("c", "c", 1400, 1000);
    l = new TLine();
    l->SetLineStyle(kDashed);
    initDrawingDone = true;
  }
  c->Divide(3, 2);
  for (int iLayer = 0; iLayer < 6; ++iLayer) {
    auto pad = c->cd(iLayer + 1);
    pad->SetRightMargin(0.15);
    hDigits[iLayer]->Draw("colz text");
    hTracklets[iLayer]->SetMarkerColor(kRed);
    hTracklets[iLayer]->SetMarkerSize(0.9);
    hTracklets[iLayer]->Draw("text same");
    hTracks[iLayer]->SetMarkerColor(kGreen);
    hTracks[iLayer]->SetMarkerSize(0.9);
    hTracks[iLayer]->Draw("text same");
    drawTrdGrid(l);
    pad->SetLogz();
  }
  c->Update();
}

bool adjustSector(o2::track::TrackParCov& t, o2::base::Propagator* prop)
{
  float alpha = geo->getAlpha();
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
    if (!prop->PropagateToXBxByBz(t, xTmp)) {
      return false;
    }
    y = t.getY();
    ++nTries;
  }
  return true;
}

int getSector(float alpha)
{
  if (alpha < 0) {
    alpha += 2.f * TMath::Pi();
  } else if (alpha >= 2.f * TMath::Pi()) {
    alpha -= 2.f * TMath::Pi();
  }
  return (int)(alpha * NSECTOR / (2.f * TMath::Pi()));
}


void trackletEfficiency(bool draw = false, int tf = -1, int trigger = -1)
{
  if (tf >= 0) {
    lastEntry = tf;
  }
  if (trigger > 0) {
    lastTrigger = trigger - 1;
  }
  if (!initDone) {
    if (!prepareInput()) {
      return;
    }
    createTrdPadHistsPerLayer();
  } else {
    for (int iLayer = 0; iLayer < 6; ++iLayer) {
      hTracklets[iLayer]->Reset();
      hDigits[iLayer]->Reset();
      hTracks[iLayer]->Reset();
      if (initDrawingDone) {
        c->Clear();
      }
    }
  }
  o2::utils::TreeStreamRedirector tStream("trdAnalysisResults.root", "recreate");
  auto prop = o2::base::Propagator::Instance();

  int countInterestingTriggers = 0;
  int countInterestingTracks = 0;
  int nTrackPoints = 0;
  int nDigitMatches = 0;
  int nTrackletMatches = 0;

  auto& ccdbmgr = o2::ccdb::BasicCCDBManager::instance();
  // auto trdHcStatus = ccdbmgr.get<o2::trd::HalfChamberStatusQC>("TRD/Calib/HalfChamberStatusQC");
  //trdHcStatus->print();

  // loop over TFs
  for (int iEntry = lastEntry; iEntry < treeDigit->GetEntries(); ++iEntry) {
    treeDigit->GetEntry(iEntry);
    treeTrklt->GetEntry(iEntry);
    treeTrks->GetEntry(iEntry);
    o2::dataformats::TFIDInfo tfid;
    if (tfids) {
      tfid = tfids->at(iEntry);
    }

    // loop over TRD triggers
    for (size_t iTrig = lastTrigger + 1; iTrig < trigIn.size(); ++iTrig) {
      const auto& trig = trigIn[iTrig];
      if (trig.getNumberOfDigits() == 0 || trig.getNumberOfTracklets() == 0) {
        continue;
      }
      float trdTriggerTimeUS = tfid.isDummy() ? trig.getBCData().bc2ns() * 1e-3 : trig.getBCData().differenceInBCMS(o2::InteractionRecord{0, tfid.firstTForbit});
      //printf("Found digits and tracklets in TF %i and trigger %lu at time %.2f us\n", iEntry, iTrig, trdTriggerTimeUS);
      //printf("In total there are %lu ITS-TPC tracks loaded for this TF\n", tracks.size());
      bool foundTrackForTrigger = false;

      // let's see if we can match tracks to digits and tracklets
      std::set<std::pair<std::tuple<int, int, int>, int>> trackMap; // filled with detector, row, col for each track point and the index of the ITS-TPC track
      std::multimap<std::tuple<int, int, int>, int> trackletMap; // filled with detector, row, col and with the index of each tracklet
      std::multimap<std::tuple<int, int, int>, int> digitMap; // filled with detector, row, col and with the index of each tracklet


      // event statistics
      int nTracks = 0;
      int nPointsTrigger = 0;

      // ITS-TPC track loop
      for (int iTrack = 0; iTrack < static_cast<int>(tracks.size()); ++iTrack) {
        const auto& trkITSTPC = tracks[iTrack];

        //printf("ITS-TPC track time: %f\n", trkITSTPC.getTimeMUS().getTimeStamp());
        if (abs(trkITSTPC.getTimeMUS().getTimeStamp() - trdTriggerTimeUS) > 2. ) {
          // track is from different interaction
          continue;
        }
        ++nTracks;
        ++countInterestingTracks;
        foundTrackForTrigger = true;
        if (draw) printf("Found ITS-TPC track with time %f us\n", trkITSTPC.getTimeMUS().getTimeStamp());
        const auto& paramOut = trkITSTPC.getParamOut();
        if (fabs(paramOut.getEta()) > 0.84 || paramOut.getPt() < minPtTrack) {
          // no chance to find tracklets for these tracks
          continue;
        }
        auto param = paramOut; // copy of const object
        for (int iLy = 0; iLy < 6; ++iLy) {
          if (!prop->PropagateToXBxByBz(param, geo->getTime0(iLy))) {
            if (draw) printf("Failed track propagation into layer %i\n", iLy);
            break;
          }
          if (!adjustSector(param, prop)) {
            if (draw) printf("Failed track rotation in layer %i\n", iLy);
            break;
          }
          if (draw) printf("Track has alpha of %f in layer %i. X(%f), Y(%f), Z(%f). Eta(%f), Pt(%f)\n", param.getAlpha(), iLy, param.getX(), param.getY(), param.getZ(), param.getEta(), param.getPt());
          auto sector = getSector(param.getAlpha());
          auto stack = geo->getStack(param.getZ(), iLy);
          if (stack < 0) {
            if (draw) printf("WARN: cannot determine stack for z = %f, layer = %i\n", param.getZ(), iLy);
            continue;
          }
          auto pp = geo->getPadPlane(iLy, stack);
          auto row = pp->getPadRowNumber(param.getZ());
          int rowMax = (stack == 2) ? 12 : 16;
          if (row < 0 || row >= rowMax) {
            if (draw) printf("WARN: row  = %i for z = %f\n", row, param.getZ());
            continue;
          }
          auto col = pp->getPadColNumber(param.getY());
          if (col < 0 || col >= 144) {
            if (draw) printf("WARN: col  = %i for y = %f\n", col, param.getY());
            continue;
          }
          ++nPointsTrigger;
          trackMap.insert(std::make_pair(std::make_tuple(geo->getDetector(iLy, stack, sector), row, col), iTrack));
          int rowGlb = stack < 3 ? row + stack * 16 : row + 44 + (stack - 3) * 16; // pad row within whole sector
          int colGlb = col + sector * 144; // pad column number from 0 to NSectors * 144
          hTracks[iLy]->SetBinContent(rowGlb+1, colGlb+1, 4);
        }
      } // end ITS-TPC track loop

      // TRD digit loop
      for (int iDigit = trig.getFirstDigit(); iDigit < trig.getFirstDigit() + trig.getNumberOfDigits(); ++iDigit) {
        const auto& digit = digitIn[iDigit];
        if (digit.isSharedDigit() || digit.getChannel() == 22) {
          continue;
        }
        int adcSum = 0;
        int adcMax = 0;
        // FIXME: the track performs the tracklet fit only for time bins 5..23
        // at least 8 clusters above 50 ADC counts are required for a tracklet to be found
        // this check should be added
        for (int iTb = 5; iTb < 24; ++iTb) {
          auto adc = digit.getADC()[iTb];
          if (adc > adcMax) {
            adcMax = adc;
          }
          adcSum += adc;
        }
        if (adcMax < minAdcCount) continue;
        int det = digit.getDetector();
        int layer = det % 6;
        int stack = (det % 30) / 6;
        int rowGlb = stack < 3 ? digit.getPadRow() + stack * 16 : digit.getPadRow() + 44 + (stack - 3) * 16; // pad row within whole sector
        int sec = det / 30;
        //if (sec == 2) std::cout << digit << std::endl;
        int colGlb = digit.getPadCol() + sec * 144; // pad column number from 0 to NSectors * 144
        digitMap.insert(std::make_pair(std::make_tuple(det, digit.getPadRow(), digit.getPadCol()), iDigit));
        hDigits[layer]->SetBinContent(rowGlb+1, colGlb+1, adcMax);
      } // end digit loop

      // TRD tracklet loop
      for (int iTrklt = trig.getFirstTracklet(); iTrklt < trig.getFirstTracklet() + trig.getNumberOfTracklets(); ++iTrklt) {
        const auto& trklt = tracklets[iTrklt];
        int det = trklt.getDetector();
        int layer = det % 6;
        int stack = (det % 30) / 6;
        int rowGlb = stack < 3 ? trklt.getPadRow() + stack * 16 : trklt.getPadRow() + 44 + (stack - 3) * 16; // pad row within whole sector
        int sec = det / 30;
        //trklt.print();
        int padLocal;
        if (doXOR) {
          padLocal = trklt.getPosition() ^ 0x80;
          if (padLocal & (1 << (NBITSTRKLPOS - 1))) {
            padLocal = -((~(padLocal - 1)) & ((1 << NBITSTRKLPOS) - 1));
          } else {
            padLocal = padLocal & ((1 << constants::NBITSTRKLPOS) - 1);
          }
        } else {
          padLocal = trklt.getPositionBinSigned();
        }
        int mcmCol = (trklt.getMCM() % constants::NMCMROBINCOL) + constants::NMCMROBINCOL * (trklt.getROB() % 2);
        float colGlb = -65.f + mcmCol * ((float) constants::NCOLMCM) + padLocal * constants::GRANULARITYTRKLPOS + 144.f * sec + 72.f;
        int colInChamber = static_cast<int>(std::round(colGlb - 144.f * sec));
        trackletMap.insert(std::make_pair(std::make_tuple(det, trklt.getPadRow(), colInChamber), iTrklt));
        hTracklets[layer]->SetBinContent(rowGlb + 1, hTracklets[layer]->GetYaxis()->FindBin(colGlb), 1);
      } // end tracklet loop

      if (foundTrackForTrigger) {
        printf("Yeah, in TF %i and trigger idx %zu we found %i digits, %i tracklets and %i track(s) making %i points at time %.2f us\n", iEntry, iTrig, trig.getNumberOfDigits(), trig.getNumberOfTracklets(), nTracks, nPointsTrigger, trdTriggerTimeUS);
        ++countInterestingTriggers;

        // search for matching digits + tracklets
        for (const auto& trackHit : trackMap) {
          const auto& trk = tracks[trackHit.second];
          // output needed only for the stream --->
          o2::dataformats::TrackTPCITS trkOut = trk;
          std::vector<Digit> digitsOut;
          std::vector<Tracklet64> trackletsOut;
          bool isMasked = false;
          // <---
          ++nTrackPoints;
          bool matchingDigit = false;
          bool matchingTracklet = false;
          auto det = std::get<0>(trackHit.first);
          auto row = std::get<1>(trackHit.first);
          auto col = std::get<2>(trackHit.first);
          if (draw) printf("\n===\nChecking track point in %02i_%i_%i, row(%i), col(%i). Seeding pT(%f), eta(%f)\n", geo->getSector(det), geo->getStack(det), geo->getLayer(det), row, col, trk.getPt(), trk.getEta());
          int halfChamber = (col >= NCOLUMN / 2) ? det * 2 + 1 : det * 2;
          // if (trdHcStatus->isMasked(halfChamber)) {
          //   isMasked = true;
          //   if (draw) printf("Track ended up in masked half chamber...\n");
          // }
          int rowMaxChamber = geo->getStack(det) == 2 ? NROWC0 : NROWC1;
          int rowMin = row - searchRoadRow;
          if (rowMin < 0) {
            rowMin = 0;
          }
          int rowMax = row + searchRoadRow;
          if (rowMax >= rowMaxChamber) {
            rowMax = rowMaxChamber - 1;
          }
          int colMin = col - searchRoadCol;
          if (colMin < 0) {
            colMin = 0;
          }
          int colMax = col + searchRoadCol;
          if (colMax >= NCOLUMN) {
            colMax = NCOLUMN - 1;
          }
          for (int iRow = rowMin; iRow <= rowMax; ++iRow) {
            for (int iCol = colMin; iCol <= colMax; ++iCol) {
              auto mapKey = std::make_tuple(det, iRow, iCol);
              auto digitRange = digitMap.equal_range(mapKey);
              std::for_each(digitRange.first, digitRange.second, [&](std::multimap<std::tuple<int, int, int>, int>::value_type& idx) {
                matchingDigit = true;
                const auto& digit = digitIn[idx.second];
                digitsOut.push_back(digit);
                if (draw) {
                  printf("Matching digit in %02i_%i_%i, row(%i), col(%i) with ADCs: ", geo->getSector(digit.getDetector()), geo->getStack(digit.getDetector()), geo->getLayer(digit.getDetector()), digit.getPadRow(), digit.getPadCol());
                  for (int iTb = 0; iTb < TIMEBINS; ++iTb) {
                    printf("[%3i]", digit.getADC()[iTb]);
                  }
                  printf("\n");
                }
              });
              auto trackletRange = trackletMap.equal_range(mapKey);
              std::for_each(trackletRange.first, trackletRange.second, [&](std::multimap<std::tuple<int, int, int>, int>::value_type& idx) {
                matchingTracklet = true;
                const auto& trklt = tracklets[idx.second];
                trackletsOut.push_back(trklt);
                if (draw)printf("Matching tracklet in %02i_%i_%i, row(%i), col(%i) with Q0(%i), Q1(%i), Q2(%i)\n", geo->getSector(trklt.getDetector()), geo->getStack(trklt.getDetector()), geo->getLayer(trklt.getDetector()), trklt.getPadRow(), iCol, trklt.getQ0(), trklt.getQ1(), trklt.getQ2());
              });
            }
          }
          if (matchingDigit) ++nDigitMatches;
          if (matchingTracklet) ++nTrackletMatches;
          tStream << "tree" <<
            "track=" << &trkOut <<
            "digits=" << &digitsOut <<
            "tracklets=" << &trackletsOut <<
            "isMasked=" << isMasked <<
            "matchedDigit=" << matchingDigit <<
            "matchedTracklet=" << matchingTracklet <<
            "\n";
        }
      }

      if (draw) {
        lastEntry = iEntry;
        lastTrigger = iTrig;
        drawResults();
        initDone = true;
        return;
      }
      initDone = true;
    } // end trigger loop
    lastTrigger = -1;
  } // end tree loop

  tStream.Close();

  printf("Found %i triggers with tracklets, digits and at least one ITS-TPC track. In total %i ITS-TPC tracks are contained in these triggers.\n", countInterestingTriggers, countInterestingTracks);
  printf("For %i track points inside the TRD we found %i matching digits and %i matching tracklets\n", nTrackPoints, nDigitMatches, nTrackletMatches);
}
