/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FXD_WIDTH_WRTR
#define FXD_WIDTH_WRTR

struct FxdWidthWrtr : Writer {

  FxdWidthWrtr(std::string const& nm, std::string const& als,
               const bool crit, const bool note) :
    Writer(nm, als, crit, note)
  {
  }

  FxdWidthWrtr & Put(const unsigned int strtIdx,
                     const unsigned int len,
                     const unsigned int ui,
                     const char algn = 'l')
  {
  }
  // int, uint, float, double, long, string, char *, bool
  //Put(...)
  // For flt/dbl, max width = 7; otherwise use full width
  // Used fixed but have PutSN (sci. not.) function that does non-fixed for
  // flt.
  // Align left by default, but accept right alignment if necessary.
}

#endif /*FXD_WIDTH_WRTR*/
