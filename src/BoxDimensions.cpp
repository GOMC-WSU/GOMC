/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "BoxDimensions.h"
#include "GeomLib.h"
#include "MoveConst.h" //For cutoff-related fail condition

using namespace geom;

void BoxDimensions::Init(config_setup::RestartSettings const& restart,
                         config_setup::Volume const& confVolume,
                         pdb_setup::Cryst1 const& cryst,
                         Forcefield const &ff)
{
  for (uint b = 0; b < BOX_TOTAL; b++) {
    rCut[b] = std::max(ff.rCut, ff.rCutCoulomb[b]);
    rCutSq[b] = rCut[b] * rCut[b];
    minVol[b] = 8.0 * rCutSq[b] * rCut[b] + 0.001;
    if(restart.enable && cryst.hasVolume) {
      axis = cryst.axis;
      real alpha = cos(cryst.cellAngle[b][0] * M_PI / 180.0);
      real beta  = cos(cryst.cellAngle[b][1] * M_PI / 180.0);
      real gamma = cos(cryst.cellAngle[b][2] * M_PI / 180.0);
      if(float(cryst.cellAngle[b][0]) == 90.0)
        alpha = 0.0;
      if(float(cryst.cellAngle[b][1]) == 90.0)
        beta = 0.0;
      if(float(cryst.cellAngle[b][2]) == 90.0)
        gamma = 0.0;
      real cosASq = alpha * alpha;
      real cosBSq = beta * beta;
      real cosGSq = gamma * gamma;
      real temp = (alpha - beta * gamma) / (sqrt(1.0 - cosGSq));
      cellBasis[b].Set(0, 1.0, 0.0, 0.0);
      cellBasis[b].Set(1, gamma, sqrt(1.0 - cosGSq), 0.0);
      cellBasis[b].Set(2, beta, temp, sqrt(1.0 - cosBSq - temp * temp));
      cellBasis[b].Scale(0, axis.Get(b).x);
      cellBasis[b].Scale(1, axis.Get(b).y);
      cellBasis[b].Scale(2, axis.Get(b).z);
    } else if(confVolume.hasVolume) {
      confVolume.axis[b].CopyRange(cellBasis[b], 0, 0, 3);
    } else {
      fprintf(stderr,
              "Error: Cell Basis not specified in PDB or in.conf files.\n");
      exit(EXIT_FAILURE);
    }

    //Print Box dimensio info
    printf("%s %-d: %-26s %6.3f %7.3f %7.3f \n",
           "Info: Box ", b, " Periodic Cell Basis 1",
           cellBasis[b].Get(0).x, cellBasis[b].Get(0).y,
           cellBasis[b].Get(0).z);
    printf("%s %-d: %-26s %6.3f %7.3f %7.3f \n",
           "Info: Box ", b, " Periodic Cell Basis 2",
           cellBasis[b].Get(1).x, cellBasis[b].Get(1).y,
           cellBasis[b].Get(1).z);
    printf("%s %-d: %-26s %6.3f %7.3f %7.3f \n\n",
           "Info: Box ", b, " Periodic Cell Basis 3",
           cellBasis[b].Get(2).x, cellBasis[b].Get(2).y,
           cellBasis[b].Get(2).z);


    axis.Set(b, cellBasis[b].Length(0), cellBasis[b].Length(1),
             cellBasis[b].Length(2));

    if(axis.Get(b).Min() < 2.0 * rCut[b]) {
      printf("Error: Cutoff value is large than half of minimum BOX%d length!\n", b);
      exit(EXIT_FAILURE);
    }
    //Find Cosine Angle of alpha, beta and gamma
    cosAngle[b][0] = Dot(cellBasis[b].Get(1), cellBasis[b].Get(2)) /
                     (axis.Get(b).y * axis.Get(b).z);
    cosAngle[b][1] = Dot(cellBasis[b].Get(0), cellBasis[b].Get(2)) /
                     (axis.Get(b).x * axis.Get(b).z);
    cosAngle[b][2] = Dot(cellBasis[b].Get(0), cellBasis[b].Get(1)) /
                     (axis.Get(b).x * axis.Get(b).y);

    volume[b] = axis.x[b] * axis.y[b] * axis.z[b];
    volInv[b] = 1.0 / volume[b];

    //normalizing unitcell
    for(uint i = 0; i < 3; i++) {
      cellBasis[b].Set(i, cellBasis[b].Get(i).Normalize());
    }
  }

  axis.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  halfAx.ScaleRange(0, BOX_TOTAL, 0.5);

  for (uint b = 0; b < BOX_TOTAL; b++) {
    //check to see if initial box size is cubic or not
    cubic[b] = ((axis.x[b] == axis.y[b]) && (axis.y[b] == axis.z[b]));
    //check to see if box is orthogonal or not
    orthogonal[b] = ((cosAngle[b][0] == 0.0) &&
                     (cosAngle[b][1] == 0.0) &&
                     (cosAngle[b][2] == 0.0));
  }
  constArea = confVolume.cstArea;
}

uint BoxDimensions::ShiftVolume(BoxDimensions & newDim, XYZ & scale,
                                const uint b, const real delta) const
{
  uint rejectState = mv::fail_state::NO_FAIL;
  real newVolume = volume[b] + delta;
  newDim.SetVolume(b, newVolume);

  //If move would shrink any box axis to be less than 2 * rCut[b], then
  //automatically reject to prevent errors.
  if ((newDim.halfAx.x[b] < rCut[b] || newDim.halfAx.y[b] < rCut[b] ||
       newDim.halfAx.z[b] < rCut[b] || newVolume < minVol[b])) {
    std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!\n";
    std::cout << "AxisDimensions: " << newDim.GetAxis(b) << std::endl;
    std::cout << "Exiting!\n";
    exit(EXIT_FAILURE);
  }
  scale = newDim.axis.Get(b) / axis.Get(b);

  return rejectState;
}

uint BoxDimensions::ExchangeVolume(BoxDimensions & newDim, XYZ * scale,
                                   const real transfer, const uint *box) const
{
  uint state = mv::fail_state::NO_FAIL;
  real vTot = GetTotVolume(box[0], box[1]);

  newDim.SetVolume(box[0], volume[box[0]] + transfer);
  newDim.SetVolume(box[1], vTot - newDim.volume[box[0]]);

  //If move would shrink any box axis to be less than 2 * rcut, then
  //automatically reject to prevent errors.
  for (uint i = 0; i < 2; i++) {
    uint b = box[i];
    scale[b] = newDim.axis.Get(b) / axis.Get(b);
    if ((newDim.halfAx.x[b] < rCut[b] || newDim.halfAx.y[b] < rCut[b] ||
         newDim.halfAx.z[b] < rCut[b] || newDim.volume[b] < minVol[b])) {
      std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!\n";
      std::cout << "AxisDimensions: " << newDim.GetAxis(b) << std::endl;
      std::cout << "Exiting!\n";
      exit(EXIT_FAILURE);
    }
  }
  return state;
}

BoxDimensions::BoxDimensions(BoxDimensions const& other) :
  axis(other.axis), halfAx(other.halfAx)
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    cellBasis[b] = XYZArray(3);
    other.cellBasis[b].CopyRange(cellBasis[b], 0, 0, 3);
    volume[b] = other.volume[b];
    volInv[b] = other.volInv[b];
    rCut[b] = other.rCut[b];
    rCutSq[b] = other.rCutSq[b];
    cubic[b] = other.cubic[b];
    orthogonal[b] = other.orthogonal[b];
    for(uint i = 0; i < 3; i++) {
      cosAngle[b][i] = other.cosAngle[b][i];
    }
  }
  constArea = other.constArea;
}

BoxDimensions& BoxDimensions::operator=(BoxDimensions const& other)
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    other.cellBasis[b].CopyRange(cellBasis[b], 0, 0, 3);
    volume[b] = other.volume[b];
    volInv[b] = other.volInv[b];
    rCut[b] = other.rCut[b];
    rCutSq[b] = other.rCutSq[b];
    cubic[b] = other.cubic[b];
    orthogonal[b] = other.orthogonal[b];
    for(uint i = 0; i < 3; i++) {
      cosAngle[b][i] = other.cosAngle[b][i];
    }
  }
  other.axis.CopyRange(axis, 0, 0, BOX_TOTAL);
  other.halfAx.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  constArea = other.constArea;
  return *this;
}

real BoxDimensions::GetTotVolume(const uint b1, const uint b2) const
{
  return (volume[b1] + volume[b2]);
}

void BoxDimensions::SetVolume(const uint b, const real vol)
{
  if(constArea) {
    real ratio = vol / volume[b];
    axis.Scale(b, 1.0, 1.0, ratio);
    halfAx.Scale(b, 1.0, 1.0, ratio);
  } else {
    real ratio = pow(vol / volume[b], (1.0 / 3.0));
    axis.Scale(b, ratio);
    halfAx.Scale(b, ratio);
  }
  volume[b] = vol;
  volInv[b] = 1.0 / vol;
}


void BoxDimensions::WrapPBC(real &x, real &y, real &z, const uint b) const
{
  WrapPBC(x, axis.x[b]);
  WrapPBC(y, axis.y[b]);
  WrapPBC(z, axis.z[b]);
}

void BoxDimensions::UnwrapPBC(real & x, real & y, real & z,
                              const uint b, XYZ const& ref) const
{
  UnwrapPBC(x, ref.x, axis.x[b], halfAx.x[b]);
  UnwrapPBC(y, ref.y, axis.y[b], halfAx.y[b]);
  UnwrapPBC(z, ref.z, axis.z[b], halfAx.z[b]);
}


//
// Note, here we can't do the fabs trick as we need to know which end
// to wrap on.
//
real BoxDimensions::WrapPBC(real& v, const real ax) const
{
  //assert(v < 2*ax);

  //if ( v > ax ) //if +, wrap out to low end
  //   v -= ax;
  //else if ( v < 0 ) //if -, wrap to high end
  //   v += ax;
  //
  // Inspired by :
  // http://graphics.stanford.edu/~seander/bithacks.html#ConditionalNegate
  //
#ifdef NO_BRANCHING_WRAP
  // Inspired by :
  // http://graphics.stanford.edu/~seander/bithacks.html#ConditionalNegate
  //
  // Advantages:
  // No branching
  //
  // Disadvantages:
  // Sometimes a couple of extra ops or a couple of extra compares.
  if (
    bool negate = (v > ax);
    real vNeg = v + (ax ^ -negate) + negate;
    return (fabs(v - halfAx) > halfAx) ? v : vNeg;
#else
  //Note: testing shows that it's most efficient to negate if true.
  //Source:
  // http://jacksondunstan.com/articles/2052
  if ( v >= ax ) //if +, wrap out to low end, on boundry will wrap to zero
    v -= ax;
  else if ( v < 0 ) //if -, wrap to high end
    v += ax;
  return v;
#endif
}

real BoxDimensions::UnwrapPBC(real& v, const real ref, const real ax,
                                const real halfAx) const
{
  //If absolute value of X dist btwn pt and ref is > 0.5 * box_axis
  //If ref > 0.5 * box_axis, add box_axis to pt (to wrap out + side)
  //If ref < 0.5 * box_axis, subtract box_axis (to wrap out - side)
  // uses bit hack to avoid branching for conditional
#ifdef NO_BRANCHING_UNWRAP
  bool negate = ( ref > halfAx );
  real vDiff = v + (ax ^ -negate) + negate;
  return (fabs(ref - v) > halfAx ) ? v : vDiff;
#else
  if (fabs(ref - v) > halfAx ) {
    //Note: testing shows that it's most efficient to negate if true.
    //Source:
    // http://jacksondunstan.com/articles/2052
    if ( ref < halfAx )
      v -= ax;
    else
      v += ax;
  }
  return v;
#endif
}

XYZ BoxDimensions::MinImage(XYZ rawVec, const uint b) const
{
  rawVec.x = MinImageSigned(rawVec.x, axis.x[b], halfAx.x[b]);
  rawVec.y = MinImageSigned(rawVec.y, axis.y[b], halfAx.y[b]);
  rawVec.z = MinImageSigned(rawVec.z, axis.z[b], halfAx.z[b]);
  return rawVec;
}

//Dist. btwn two points, accounting for PBC, on an individual axis
//
//Throws out sign (as per Brock's suggestion) as we don't care about it
//and thus can eliminate a branch and (potentially) one compare.
//
real BoxDimensions::MinImage
(real& raw, const real ax, const real halfAx) const
{
  raw = fabs(raw);
  //If shorter over periodic boundary, get that dist.
#ifdef NO_BRANCHING_MIN_IMAGE
  rawDiff = ax - raw;
  return (raw > halfAx) ? rawDiff : raw;
#else
  if (raw > halfAx)
    raw = ax - raw;
  return raw; //...just pass back if distance is already minimum
#endif
}

real BoxDimensions::MinImageSigned(real raw, real ax, real halfAx) const
{
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}
