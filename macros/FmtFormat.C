

#ifndef TRD_PAYGROUND_FMT_FORMAT_C
#define TRD_PAYGROUND_FMT_FORMAT_C

#include <fmt/format.h>
#include <map>

#include <DataFormatsTRD/Digit.h>
#include <DataFormatsTRD/Tracklet64.h>

template <>
struct fmt::formatter<o2::trd::Digit>
{
  // Presentation format: 'c' - coordinate only, 'a' - all (not yet implemented).
  char presentation = 'c';

  // Parses format specifications of the form ['f' | 'e'].
  constexpr auto
  parse(format_parse_context &ctx) -> decltype(ctx.begin())
  {
    // [ctx.begin(), ctx.end()) is a character range that contains a part of
    // the format string starting from the format specifications to be parsed,
    // e.g. in
    //
    //   fmt::format("{:f} - point of interest", point{1, 2});
    //
    // the range will contain "f} - point of interest". The formatter should
    // parse specifiers until '}' or the end of the range. In this example
    // the formatter should parse the 'f' specifier and return an iterator
    // pointing to '}'.

    // Parse the presentation format and store it in the formatter:
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && (*it == 'c' || *it == 'r' || *it == 'm'))
      presentation = *it++;

    // Check if reached the end of the range:
    if (it != end && *it != '}')
      throw format_error("invalid format");

    // Return an iterator past the end of the parsed range:
    return it;
  }

  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  auto format(const o2::trd::Digit &digit, FormatContext &ctx) -> decltype(ctx.out())
  {

    int sector = digit.getDetector() / o2::trd::constants::NCHAMBERPERSEC;
    int stack = (digit.getDetector() % o2::trd::constants::NCHAMBERPERSEC) / o2::trd::constants::NLAYER;
    int layer = digit.getDetector() % o2::trd::constants::NLAYER;

                // ctx.out() is an output iterator to write to.
                switch (presentation)
    {
    case 'c':
      return format_to(ctx.out(), "{:02d}_{:d}_{:d} {:d}:{:02d}.{:02d} ( {:2d} / {:3d} )",
                       sector, stack, layer, 
                       digit.getROB(), digit.getMCM(), digit.getChannel(), 
                       digit.getPadRow(), digit.getPadCol());

    case 'r':
      return format_to(ctx.out(), "{:02d}_{:d}_{:d} row {:d}",
                       sector, stack, layer, digit.getPadRow());

    case 'm':
      return format_to(ctx.out(), "{:02d}_{:d}_{:d} {:d}:{:02d}",
                       sector, stack, layer, digit.getROB(), digit.getMCM());

    default:
      return format_to(ctx.out(), "{:02d}_{:d}_{:d} {:d}:{:02d}",
                       sector, stack, layer, digit.getROB(), digit.getMCM());
    }
  }
};


template <>
struct fmt::formatter<o2::trd::Tracklet64>
{
  // Presentation format: 'c' - coordinate only, 'a' - all (not yet implemented).
  char presentation = 'c';

  // Parses format specifications of the form ['f' | 'e'].
  constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin())
  {
    // [ctx.begin(), ctx.end()) is a character range that contains a part of
    // the format string starting from the format specifications to be parsed,
    // e.g. in
    //
    //   fmt::format("{:f} - point of interest", point{1, 2});
    //
    // the range will contain "f} - point of interest". The formatter should
    // parse specifiers until '}' or the end of the range. In this example
    // the formatter should parse the 'f' specifier and return an iterator
    // pointing to '}'.

    // Parse the presentation format and store it in the formatter:
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && (*it == 'c' || *it == 'a'))
      presentation = *it++;

    // Check if reached the end of the range:
    if (it != end && *it != '}')
      throw format_error("invalid format");

    // Return an iterator past the end of the parsed range:
    return it;
  }

  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  auto format(const o2::trd::Tracklet64 &tracklet, FormatContext &ctx) -> decltype(ctx.out())
  {

    int sector = tracklet.getDetector() / o2::trd::constants::NCHAMBERPERSEC;
    int stack = (tracklet.getDetector() % o2::trd::constants::NCHAMBERPERSEC) / o2::trd::constants::NLAYER;
    int layer = tracklet.getDetector() % o2::trd::constants::NLAYER;

    // ctx.out() is an output iterator to write to.
    switch (presentation)
    {
    case 'c':
      return format_to(ctx.out(), "{:02d}_{:d}_{:d} {:d}:{:02d} pos={:d} slope={:d}",
                       sector, stack, layer,
                       tracklet.getROB(), tracklet.getMCM(), 
                       tracklet.getPosition(), tracklet.getSlope());

    default:
      return format_to(ctx.out(), "ASD"); //"{:02d}_{d}_{d} {d}:{:02d}",
                      //  sector, stack, layer, tracklet.getROB(), tracklet.getMCM());
    }
  }
};

#endif
