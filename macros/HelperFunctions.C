
#ifndef TRD_PLAYGROUND_HELPER_FUNCTIONS_H
#define TRD_PLAYGROUND_HELPER_FUNCTIONS_H

#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/HelperMethods.h"
using namespace o2::trd::constants;

int getMCMCol(int irob, int imcm)
{
  return (imcm % NMCMROBINCOL) + NMCMROBINCOL * (irob % 2);
}

float PadPositionMCM(o2::trd::Tracklet64 &tracklet)
{
  return 12.0 - (tracklet.getPositionBinSigned() * GRANULARITYTRKLPOS);
}

float PadPosition(o2::trd::Tracklet64 &tracklet)
{
  // int padLocalBin = tracklet.getPosition() ^ 0x80;
  float padMCM = PadPositionMCM(tracklet);
  int mcmCol = getMCMCol(tracklet.getROB(), tracklet.getMCM());
  return float( (mcmCol+1) * NCOLMCM ) + 2.0 - padMCM;
}

float UncalibratedPad(o2::trd::Tracklet64 &tracklet)
{
  float y = tracklet.getUncalibratedY();
  int mcmCol = (tracklet.getMCM() % NMCMROBINCOL) + NMCMROBINCOL * (tracklet.getROB() % 2);
  // one pad column has 144 pads, the offset of -63 is the center of the first MCM in that column
  // which is connected to the pads -63 - 9 = -72 to -63 + 9 = -54
  // float offset = -63.f + ((float)NCOLMCM) * mcmCol;
  float padWidth = 0.635f + 0.03f * (tracklet.getDetector() % NLAYER);
  return y / padWidth + 71.0;
}

float Slope(o2::trd::Tracklet64 &trkl)
{
  return trkl.getSlopeBinSigned() * GRANULARITYTRKLSLOPE / ADDBITSHIFTSLOPE;
}

#endif