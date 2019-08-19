/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
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
      axis[b][0] = cryst.axis.x[b];
      axis[b][1] = cryst.axis.y[b];
      axis[b][2] = cryst.axis.z[b];
      double alpha = cos(cryst.cellAngle[b][0] * M_PI / 180.0);
      double beta  = cos(cryst.cellAngle[b][1] * M_PI / 180.0);
      double gamma = cos(cryst.cellAngle[b][2] * M_PI / 180.0);
      if(float(cryst.cellAngle[b][0]) == 90.0)
        alpha = 0.0;
      if(float(cryst.cellAngle[b][1]) == 90.0)
        beta = 0.0;
      if(float(cryst.cellAngle[b][2]) == 90.0)
        gamma = 0.0;
      double cosASq = alpha * alpha;
      double cosBSq = beta * beta;
      double cosGSq = gamma * gamma;
      double temp = (alpha - beta * gamma) / (sqrt(1.0 - cosGSq));
      cellBasis[b].Set(0, 1.0, 0.0, 0.0);
      cellBasis[b].Set(1, gamma, sqrt(1.0 - cosGSq), 0.0);
      cellBasis[b].Set(2, beta, temp, sqrt(1.0 - cosBSq - temp * temp));
      cellBasis[b].Scale(0, axis[b][0]);
      cellBasis[b].Scale(1, axis[b][1]);
      cellBasis[b].Scale(2, axis[b][2]);
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


    axis[b][0] = cellBasis[b].Length(0);
    axis[b][1] = cellBasis[b].Length(1);
    axis[b][2] = cellBasis[b].Length(2);

    if(std::min(std::min(axis[b][0], axis[b][1]), axis[b][2]) < 2.0 * rCut[b]) {
      printf("Error: Cutoff value is large than half of minimum BOX%d length!\n", b);
      exit(EXIT_FAILURE);
    }
    //Find Cosine Angle of alpha, beta and gamma
    cosAngle[b][0] = Dot(cellBasis[b].Get(1), cellBasis[b].Get(2)) /
                     (axis[b][1] * axis[b][2]);
    cosAngle[b][1] = Dot(cellBasis[b].Get(0), cellBasis[b].Get(2)) /
                     (axis[b][0] * axis[b][2]);
    cosAngle[b][2] = Dot(cellBasis[b].Get(0), cellBasis[b].Get(1)) /
                     (axis[b][0] * axis[b][1]);

    volume[b] = axis[b][0] * axis[b][1] * axis[b][2];
    volInv[b] = 1.0 / volume[b];

    //normalizing unitcell
    for(uint i = 0; i < 3; i++) {
      cellBasis[b].Set(i, cellBasis[b].Get(i).Normalize());
    }
  }

  for (uint b = 0; b < BOX_TOTAL; b++) {
    halfAx[b][0] = axis[b][0] * 0.5;
    halfAx[b][1] = axis[b][1] * 0.5;
    halfAx[b][2] = axis[b][2] * 0.5;
  }

  for (uint b = 0; b < BOX_TOTAL; b++) {
    //check to see if initial box size is cubic or not
    cubic[b] = ((axis[b][0] == axis[b][1]) && (axis[b][1] == axis[b][2]));
    //check to see if box is orthogonal or not
    orthogonal[b] = ((cosAngle[b][0] == 0.0) &&
                     (cosAngle[b][1] == 0.0) &&
                     (cosAngle[b][2] == 0.0));
  }
  constArea = confVolume.cstArea;
}

uint BoxDimensions::ShiftVolume(BoxDimensions & newDim, XYZ & scale,
                                const uint b, const double delta) const
{
  uint rejectState = mv::fail_state::NO_FAIL;
  double newVolume = volume[b] + delta;
  newDim.SetVolume(b, newVolume);

  //If move would shrink any box axis to be less than 2 * rCut[b], then
  //automatically reject to prevent errors.
  if ((newDim.halfAx[b][0] < rCut[b] || newDim.halfAx[b][1] < rCut[b] ||
       newDim.halfAx[b][2] < rCut[b] || newVolume < minVol[b])) {
    std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!\n";
    std::cout << "AxisDimensions: " << newDim.GetAxis(b) << std::endl;
    std::cout << "Exiting!\n";
    exit(EXIT_FAILURE);
  }
  scale.x = newDim.axis[b][0] / axis[b][0];
  scale.y = newDim.axis[b][1] / axis[b][1];
  scale.z = newDim.axis[b][2] / axis[b][2];

  return rejectState;
}

uint BoxDimensions::ExchangeVolume(BoxDimensions & newDim, XYZ * scale,
                                  const double transfer, const uint *box) const
{
  uint state = mv::fail_state::NO_FAIL;
  double vTot = GetTotVolume(box[0], box[1]);

  newDim.SetVolume(box[0], volume[box[0]] + transfer);
  newDim.SetVolume(box[1], vTot - newDim.volume[box[0]]);

  //If move would shrink any box axis to be less than 2 * rcut, then
  //automatically reject to prevent errors.
  for (uint i = 0; i < 2; i++) {
    uint b = box[i];
    scale[b].x = newDim.axis[b][0] / axis[b][0];
    scale[b].y = newDim.axis[b][1] / axis[b][1];
    scale[b].z = newDim.axis[b][2] / axis[b][2];
    if ((newDim.halfAx[b][0] < rCut[b] || newDim.halfAx[b][1] < rCut[b] ||
         newDim.halfAx[b][2] < rCut[b] || newDim.volume[b] < minVol[b])) {
      std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!\n";
      std::cout << "AxisDimensions: " << newDim.GetAxis(b) << std::endl;
      std::cout << "Exiting!\n";
      exit(EXIT_FAILURE);
    }
  }
  return state;
}

BoxDimensions::BoxDimensions(BoxDimensions const& other)
{
  memcpy(axis, other.axis, BOX_TOTAL * 3 * sizeof(double));
  memcpy(halfAx, other.halfAx, BOX_TOTAL * 3 * sizeof(double));

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
  memcpy(axis, other.axis, BOX_TOTAL * 3 * sizeof(double));
  memcpy(halfAx, other.halfAx, BOX_TOTAL * 3 * sizeof(double));
  constArea = other.constArea;
  return *this;
}

double BoxDimensions::GetTotVolume(const uint b1, const uint b2) const
{
  return (volume[b1] + volume[b2]);
}

void BoxDimensions::SetVolume(const uint b, const double vol)
{
  if(constArea) {
    double ratio = vol / volume[b];
    axis[b][2] *= ratio;
    halfAx[b][2] *= ratio;
  } else {
    double ratio = pow(vol / volume[b], (1.0 / 3.0));
    axis[b][0] *= ratio;
    axis[b][1] *= ratio;
    axis[b][2] *= ratio;
    halfAx[b][0] *= ratio;
    halfAx[b][1] *= ratio;
    halfAx[b][2] *= ratio;
  }
  volume[b] = vol;
  volInv[b] = 1.0 / vol;
}


void BoxDimensions::WrapPBC(double &x, double &y, double &z, const uint b) const
{
  WrapPBC(x, axis[b][0]);
  WrapPBC(y, axis[b][1]);
  WrapPBC(z, axis[b][2]);
}

void BoxDimensions::UnwrapPBC(double & x, double & y, double & z,
                              const uint b, XYZ const& ref) const
{
  UnwrapPBC(x, ref.x, axis[b][0], halfAx[b][0]);
  UnwrapPBC(y, ref.y, axis[b][1], halfAx[b][1]);
  UnwrapPBC(z, ref.z, axis[b][2], halfAx[b][2]);
}


//
// Note, here we can't do the fabs trick as we need to know which end
// to wrap on.
//
double BoxDimensions::WrapPBC(double& v, const double ax) const
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
    double vNeg = v + (ax ^ -negate) + negate;
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

double BoxDimensions::UnwrapPBC(double& v, const double ref, const double ax,
                                const double halfAx) const
{
  //If absolute value of X dist btwn pt and ref is > 0.5 * box_axis
  //If ref > 0.5 * box_axis, add box_axis to pt (to wrap out + side)
  //If ref < 0.5 * box_axis, subtract box_axis (to wrap out - side)
  // uses bit hack to avoid branching for conditional
#ifdef NO_BRANCHING_UNWRAP
  bool negate = ( ref > halfAx );
  double vDiff = v + (ax ^ -negate) + negate;
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
  // It is way faster to do it here than to call MinImageSigned 3 times
  if (rawVec.x > halfAx[b][0])
    rawVec.x = rawVec.x - axis[b][0];
  else if (rawVec.x < -halfAx[b][0])
    rawVec.x = rawVec.x + axis[b][0];
  if (rawVec.y > halfAx[b][1])
    rawVec.y = rawVec.y - axis[b][1];
  else if (rawVec.y < -halfAx[b][1])
    rawVec.y = rawVec.y + axis[b][1];
  if (rawVec.z > halfAx[b][2])
    rawVec.z = rawVec.z - axis[b][2];
  else if (rawVec.z < -halfAx[b][2])
    rawVec.z = rawVec.z + axis[b][2];
  return rawVec;
}

//Dist. btwn two points, accounting for PBC, on an individual axis
//
//Throws out sign (as per Brock's suggestion) as we don't care about it
//and thus can eliminate a branch and (potentially) one compare.
//
double BoxDimensions::MinImage
(double& raw, const double ax, const double halfAx) const
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

double BoxDimensions::MinImageSigned(double raw, double ax, double halfAx) const
{
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}
