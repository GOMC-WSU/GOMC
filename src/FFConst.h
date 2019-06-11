/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_CONST_H
#define FF_CONST_H

#include <string> //For aliases

#include "BasicTypes.h" //For uint
#include "GeomLib.h" //For RadToDeg

namespace ff
{
static const real DENSITY_N_PER_V_TO_G_ML = 1.660538921;

namespace part
{
extern const std::string WILD;  // "X"
const uint nm_len = 3;
const uint lj_n = 12;
const uint lj_Cn = 4;
//Flag to indicate whether particle is normal kind or should be handled
//via Mie potential.
namespace variety
{
const uint lj = 0;
const uint mie = 1;
const uint cnt = 2;
}
}

namespace dih
{
namespace rot
{
const uint cis = 0;
const uint trans = 1;
inline uint GetRotKind(const real phi)
{
  real ang = geom::RadToDeg(phi);
  return ((ang > 90 || ang < -90) ? ff::dih::rot::trans : ff::dih::rot::cis);
}
}

namespace kind
{
const uint MULTI = 0;
const uint POWER = 1;
const uint QUAD = 2;
}
}

namespace kind
{
namespace charmm
{
const uint nb_lj = 0;
const uint bnd = 1;
const uint ang = 2;
const uint dih = 3;
const uint impr = 4;
}
namespace exotic
{
const uint dih_pwr = 0;
const uint dih_quad = 1;
const uint nb_mie = 2;
}
const uint FILE_CHARMM = 0;
const uint FILE_EXOTIC = 1;
const uint FILE_KIND_COUNT = 2;
}
}

#endif /*FF_CONST_H*/
