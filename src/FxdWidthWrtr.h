/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef FXD_WIDTH_WRTR_H
#define FXD_WIDTH_WRTR_H

struct FxdWidthWrtr : Writer {
  FxdWidthWrtr(std::string const &nm, std::string const &als, const bool crit,
               const bool note)
      : Writer(nm, als, crit, note) {}

  FxdWidthWrtr &Put(const unsigned int strtIdx, const unsigned int len,
                    const unsigned int ui, const char algn = 'l') {}
  // int, uint, float, double, long, string, char *, bool
  // Put(...)
  // For flt/dbl, max width = 7; otherwise use full width
  // Used fixed but have PutSN (sci. not.) function that does non-fixed for
  // flt.
  // Align left by default, but accept right alignment if necessary.
}

#endif /*FXD_WIDTH_WRTR_H*/
