/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.1
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "MoveConst.h" //For cutoff-related fail condition

void BoxDimensions::Init(config_setup::RestartSettings const& restart,
                         config_setup::Volume const& confVolume,
                         pdb_setup::Cryst1 const& cryst,
                         double rc, double rcSq)
{
  rCut = rc;
  rCutSq = rcSq;
  for (uint b = 0; b < BOX_TOTAL; b++)
  {
    //if (restart.enable && cryst.hasVolume)
    //  axis = cryst.axis;
    if(confVolume.hasVolume)
      confVolume.axis[b].CopyRange(cellBasis[b], 0, 0, 3);
    else
    {
      fprintf(stderr,
	      "Error: Cell Basis not specified in PDB or in.dat files.\n");
      exit(EXIT_FAILURE);
    }
    
    axis.Set(b, cellBasis[b].Length(0), cellBasis[b].Length(1),
	     cellBasis[b].Length(2));
    //Find Cosine Angle of alpha, beta and gamma
    cosAngle[b][0] = DotProduct(cellBasis[b].Get(1), cellBasis[b].Get(2)) /
      (axis.Get(b).y * axis.Get(b).z);
    cosAngle[b][1] = DotProduct(cellBasis[b].Get(0), cellBasis[b].Get(2)) /
      (axis.Get(b).x * axis.Get(b).z);
    cosAngle[b][2] = DotProduct(cellBasis[b].Get(0), cellBasis[b].Get(1)) /
      (axis.Get(b).x * axis.Get(b).y);	   

    volume[b] = axis.x[b] * axis.y[b] * axis.z[b];
    volInv[b] = 1.0 / volume[b];
  }

  axis.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  halfAx.ScaleRange(0, BOX_TOTAL, 0.5);

  for (uint b = 0; b < BOX_TOTAL; b++)
  {
    //check to see if initial box size is cubic or not
    cubic[b] = ((axis.x[b] == axis.y[b]) && (axis.y[b] == axis.z[b]));
    //check to see if box is orthogonal or not
    orthogonal[b] = ((cosAngle[b][0] == 0.0) &&
		     (cosAngle[b][1] == 0.0) &&
		     (cosAngle[b][2] == 0.0));
  }
  constArea = confVolume.cstArea;
}

uint BoxDimensions::ShiftVolume
(BoxDimensions & newDim, XYZ & scale, const uint b, const double delta) const
{
  uint rejectState = mv::fail_state::NO_FAIL;
  double newVolume = volume[b] + delta;

  newDim.SetVolume(b, newVolume);

  //If move would shrink any box axis to be less than 2 * rcut, then
  //automatically reject to prevent errors.
  if ((newDim.halfAx.x[b] < rCut || newDim.halfAx.y[b] < rCut ||
       newDim.halfAx.z[b] < rCut))
  {
    std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!"
	      << std::endl;
    rejectState = mv::fail_state::VOL_TRANS_WOULD_SHRINK_BOX_BELOW_CUTOFF;
  }
  scale = newDim.axis.Get(b) / axis.Get(b);

  return rejectState;
}

uint BoxDimensions::ExchangeVolume
(BoxDimensions & newDim, XYZ * scale, const double transfer) const
{
  uint state = mv::fail_state::NO_FAIL;
  double vTot = GetTotVolume();
  newDim = *this;

  newDim.SetVolume(0, volume[0] + transfer);
  newDim.SetVolume(1, vTot - newDim.volume[0]);

  //If move would shrink any box axis to be less than 2 * rcut, then
  //automatically reject to prevent errors.
  for (uint b = 0; b < BOX_TOTAL && state == mv::fail_state::NO_FAIL; b++)
  {
    scale[b] = newDim.axis.Get(b) / axis.Get(b);
    if ((newDim.halfAx.x[b] < rCut || newDim.halfAx.y[b] < rCut ||
	 newDim.halfAx.z[b] < rCut))
    {
      std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!"
	      << std::endl;
      state = state && mv::fail_state::VOL_TRANS_WOULD_SHRINK_BOX_BELOW_CUTOFF;
    }
  }
  return state;
}

BoxDimensions::BoxDimensions(BoxDimensions const& other) :
  axis(other.axis), halfAx(other.halfAx)
{
  for (uint b = 0; b < BOX_TOTAL; ++b)
  {
    cellBasis[b] = XYZArray(3);
    other.cellBasis[b].CopyRange(cellBasis[b], 0, 0, 3);
    volume[b] = other.volume[b];
    volInv[b] = other.volInv[b];
    cubic[b] = other.cubic[b];
    orthogonal[b] = other.orthogonal[b];
    for(uint i = 0; i < 0; i++)
    {
      cosAngle[b][i] = other.cosAngle[b][i];
    }
  }
  rCut = other.rCut;
  rCutSq = other.rCutSq;
  constArea = other.constArea;
}

BoxDimensions& BoxDimensions::operator=(BoxDimensions const& other)
{
  for (uint b = 0; b < BOX_TOTAL; ++b)
  {
    other.cellBasis[b].CopyRange(cellBasis[b], 0, 0, 3);
    volume[b] = other.volume[b];
    volInv[b] = other.volInv[b];
    cubic[b] = other.cubic[b];
    orthogonal[b] = other.orthogonal[b];
    for(uint i = 0; i < 0; i++)
    {
      cosAngle[b][i] = other.cosAngle[b][i];
    }
  }
  other.axis.CopyRange(axis, 0, 0, BOX_TOTAL);
  other.halfAx.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  rCut = other.rCut;
  rCutSq = other.rCutSq;
  constArea = other.constArea;

  return *this;
}

double BoxDimensions::GetTotVolume() const
{
  double sum = 0.0;
  for (uint b = 0; b < BOX_TOTAL; b++)
    sum += volume[b];
  return sum;
}


void BoxDimensions::SetVolume(const uint b, const double vol)
{
   if (cubic[b] && !constArea)
   {
      double newAxX_b = pow(vol, (1.0/3.0));
      XYZ newAx_b(newAxX_b, newAxX_b, newAxX_b);
      axis.Set(b, newAx_b);
      newAx_b *= 0.5;
      halfAx.Set(b, newAx_b);
   }
   else if (constArea)
   {
     double area = axis.x[b] * axis.y[b];
     double newAxX_b = vol / area;
     XYZ newAx_b(axis.x[b], axis.y[b], newAxX_b);
     axis.Set(b, newAx_b);
     newAx_b *= 0.5;
     halfAx.Set(b, newAx_b);
   }
   else
   {
     double z_x = axis.z[b] / axis.x[b];
     double y_x = axis.y[b] / axis.x[b];
     double newAxX_b = pow(vol / (z_x * y_x), (1.0/3.0));
     XYZ newAx_b(newAxX_b, y_x * newAxX_b, newAxX_b * z_x);
     axis.Set(b, newAx_b);
     newAx_b *= 0.5;
     halfAx.Set(b, newAx_b);
   }

   cellBasis[b].Set(0, XYZ(axis.x[b], 0, 0));
   cellBasis[b].Set(1, XYZ(0, axis.y[b], 0));
   cellBasis[b].Set(2, XYZ(0, 0, axis.z[b]));
   volume[b] = vol;
   volInv[b] = 1.0/vol;

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
double BoxDimensions::MinImage
(double& raw, const double ax, const double halfAx) const
{
  raw = fabs(raw);
  //If shorter over periodic boundary, get that dist.
#ifdef NO_BRANCHING_MIN_IMAGE
  rawDiff = ax-raw;
  return (raw > halfAx)?rawDiff:raw;
#else
  if (raw > halfAx)
    raw = ax-raw;
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

bool BoxDimensions::InRcut(double & distSq, XYZ & dist,
			   XYZArray const& arr, const uint i,
			   const uint j, const uint b) const
{
  dist = MinImage(arr.Difference(i, j), b);
  distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
  return (rCutSq > distSq);
}


bool BoxDimensions::InRcut(double & distSq, XYZ & dist,
				  XYZArray const& arr1, const uint i,
				  XYZArray const& arr2, const uint j,
				  const uint b) const
{
  dist = MinImage(arr1.Difference(i, arr2, j), b);
  distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
  return (rCutSq > distSq);
}

bool BoxDimensions::InRcut(double & distSq, XYZArray const& arr,
				  const uint i, const uint j,
				  const uint b) const
{
  XYZ dist = MinImage(arr.Difference(i, j), b);
  distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
  return (rCutSq > distSq);
}


bool BoxDimensions::InRcut(double & distSq, XYZArray const& arr1,
				  const uint i, XYZArray const& arr2,
				  const uint j, const uint b) const
{
  XYZ dist = MinImage(arr1.Difference(i, arr2, j), b);
  distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
  return (rCutSq > distSq);
}


void BoxDimensions::GetDistSq(double & distSq, XYZArray const& arr1,
                                     const uint i, XYZArray const& arr2,
				     const uint j, const uint b) const
{
  XYZ dist = MinImage(arr1.Difference(i, arr2, j), b);
  distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
}

void BoxDimensions::GetDistSq(double & distSq, XYZArray const& arr,
				     const uint i, const uint j,
				     const uint b) const
{
  XYZ dist = MinImage(arr.Difference(i, j), b);
  distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
}

//Wrap one coordinate.
XYZ BoxDimensions::WrapPBC(XYZ rawPos, const uint b) const
{
  WrapPBC(rawPos.x, rawPos.y, rawPos.z, b);
  return rawPos;
}

//Wrap all coordinates in object.
void BoxDimensions::WrapPBC(XYZArray & arr, const uint b) const
{
  for (uint i = 0; i < arr.count; i++)
    WrapPBC(arr.x[i], arr.y[i], arr.z[i], b);
}

//Unwrap all coordinates in object.
void BoxDimensions::UnwrapPBC(XYZArray & arr, const uint b, XYZ
                                     const& ref) const
{
  for (uint i = 0; i < arr.count; i++)
    UnwrapPBC(arr.x[i], arr.y[i], arr.z[i], b, ref);
}

//Wrap range of coordinates in object
void BoxDimensions::WrapPBC
(XYZArray & arr, const uint start, const uint stop, const uint b) const
{
  for (uint i = start; i < stop; i++)
    WrapPBC(arr.x[i], arr.y[i], arr.z[i], b);
}

void BoxDimensions::WrapPBC(double &x, double &y, double &z,
				   const uint b) const
{
  WrapPBC(x, axis.x[b]);
  WrapPBC(y, axis.y[b]);
  WrapPBC(z, axis.z[b]);
}

void BoxDimensions::UnwrapPBC(double & x, double & y, double & z,
				     const uint b, XYZ const& ref) const
{
  UnwrapPBC(x, ref.x, axis.x[b], halfAx.x[b]);
  UnwrapPBC(y, ref.y, axis.y[b], halfAx.y[b]);
  UnwrapPBC(z, ref.z, axis.z[b], halfAx.z[b]);
}

//Unwrap range of coordinates in object
void BoxDimensions::UnwrapPBC(XYZArray & arr, const uint start,
                                     const uint stop, const uint b,
                                     XYZ const& ref) const
{
  for (uint i = start; i < stop; i++)
    UnwrapPBC(arr.x[i], arr.y[i], arr.z[i], b, ref);
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
    double vNeg = v+(ax ^-negate)+negate;
    return (fabs(v-halfAx)>halfAx)?v:vNeg;
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

//Unwrap one coordinate.
XYZ BoxDimensions::UnwrapPBC(XYZ & rawPos, const uint b,
                                    XYZ const& ref) const
{
  UnwrapPBC(rawPos.x, rawPos.y, rawPos.z, b, ref);
  return rawPos;
}

double BoxDimensions::UnwrapPBC
(double& v, const double ref, const double ax, const double halfAx) const
{
  //If absolute value of X dist btwn pt and ref is > 0.5 * box_axis
  //If ref > 0.5 * box_axis, add box_axis to pt (to wrap out + side)
  //If ref < 0.5 * box_axis, subtract box_axis (to wrap out - side)
  // uses bit hack to avoid branching for conditional
#ifdef NO_BRANCHING_UNWRAP
  bool negate = ( ref > halfAx );
  double vDiff = v+(ax^-negate)+negate;
  return (fabs(ref-v) > halfAx )?v:vDiff;
#else
  if (fabs(ref-v) > halfAx )
  {
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


//Calculate dot product
double BoxDimensions::DotProduct(const uint atom, double kx,
                                        double ky, double kz,
                                        const XYZArray &Coords) const
{
  double x = Coords.x[atom], y = Coords.y[atom], z = Coords.z[atom];
  return(x * kx + y * ky + z * kz);
}

//Calculate dot product
double BoxDimensions::DotProduct(const XYZ &A,const XYZ &B) const
{
  return(A.x * B.x + A.y * B.y + A.z * B.z);
}

//Calculate transform
XYZ BoxDimensions::TransformSlant(const XYZ &A,const uint b) const
{
  return A;
}

//Calculate transform
XYZ BoxDimensions::TransformUnSlant(const XYZ &A, const uint b) const
{
  return A;
}
