/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.1
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "BoxDimensionsNonOrth.h"
#include "BoxDimensions.h"
#include "MoveConst.h" //For cutoff-related fail condition

void BoxDimensionsNonOrth::Init(config_setup::RestartSettings const& restart,
                                config_setup::Volume const& confVolume,
                                pdb_setup::Cryst1 const& cryst,
                                double rc, double rcSq)
{
  rCut = rc;
  rCutSq = rcSq;
  for (uint b = 0; b < BOX_TOTAL; b++) {
    if(restart.enable && cryst.hasVolume) {
      axis = cryst.axis;
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
      cellBasis[b].Scale(0, axis.Get(b).x);
      cellBasis[b].Scale(1, axis.Get(b).y);
      cellBasis[b].Scale(2, axis.Get(b).z);
    } else if(confVolume.hasVolume) {
      confVolume.axis[b].CopyRange(cellBasis[b], 0, 0, 3);
    } else {
      fprintf(stderr,
              "Error: Cell Basis not specified in PDB or in.dat files.\n");
      exit(EXIT_FAILURE);
    }
    //Find the length of a, b, c
    cellLength.Set(b, cellBasis[b].Length(0), cellBasis[b].Length(1),
                   cellBasis[b].Length(2));
    //Find Cosine Angle of alpha, beta and gamma
    cosAngle[b][0] = DotProduct(cellBasis[b].Get(1), cellBasis[b].Get(2)) /
                     (cellLength.Get(b).y * cellLength.Get(b).z);
    cosAngle[b][1] = DotProduct(cellBasis[b].Get(0), cellBasis[b].Get(2)) /
                     (cellLength.Get(b).x * cellLength.Get(b).z);
    cosAngle[b][2] = DotProduct(cellBasis[b].Get(0), cellBasis[b].Get(1)) /
                     (cellLength.Get(b).x * cellLength.Get(b).y);
    //Calculate Cross Product
    XYZ axb = CrossProduct(cellBasis[b].Get(0), cellBasis[b].Get(1));
    XYZ bxc = CrossProduct(cellBasis[b].Get(1), cellBasis[b].Get(2));
    XYZ cxa = CrossProduct(cellBasis[b].Get(2), cellBasis[b].Get(0));
    //Calculate volume = A.(B x C)
    volume[b] = abs(DotProduct(cellBasis[b].Get(0), bxc));
    volInv[b] = 1.0 / volume[b];
    //normalizing unitcell
    for(uint i = 0; i < 3; i++) {
      cellBasis[b].Set(i, cellBasis[b].Get(i).Normalize());
    }
    //Calculate the adjoint and determinant
    double det = cellBasis[b].AdjointMatrix(cellBasis_Inv[b]);
    //Calculate the inverse matrix of cell basis
    cellBasis_Inv[b].ScaleRange(0, 3, 1.0 / det);
    //Set the axis with unslant cell box
    //XYZ unslant = TransformUnSlant(cellLength.Get(b), b);
    //axis.Set(b, unslant.x, unslant.y, unslant.z);
    axis.Set(b, cellLength[b]);
  }
  //Set half axis
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

void BoxDimensionsNonOrth::CalcCellDimensions(const uint b)
{
  //normalizing unitcell
  for(uint i = 0; i < 3; i++) {
    cellBasis[b].Set(i, cellBasis[b].Get(i).Normalize());
  }
  //Calculate the adjoint and determinant
  double det = cellBasis[b].AdjointMatrix(cellBasis_Inv[b]);
  //Calculate the inverse matrix of cell basis
  cellBasis_Inv[b].ScaleRange(0, 3, 1.0 / det);
}


BoxDimensionsNonOrth& BoxDimensionsNonOrth::operator=(BoxDimensionsNonOrth const& other)
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    other.cellBasis[b].CopyRange(cellBasis[b], 0, 0, 3);
    other.cellBasis_Inv[b].CopyRange(cellBasis_Inv[b], 0, 0, 3);
    volume[b] = other.volume[b];
    volInv[b] = other.volInv[b];
    cubic[b] = other.cubic[b];
    orthogonal[b] = other.orthogonal[b];
    for(uint i = 0; i < 0; i++) {
      cosAngle[b][i] = other.cosAngle[b][i];
    }
  }
  other.axis.CopyRange(axis, 0, 0, BOX_TOTAL);
  other.halfAx.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  other.cellLength.CopyRange(cellLength, 0, 0, BOX_TOTAL);
  rCut = other.rCut;
  rCutSq = other.rCutSq;
  constArea = other.constArea;
  return *this;
}


uint BoxDimensionsNonOrth::ShiftVolume(BoxDimensionsNonOrth & newDim,
                                       XYZ & scale, const uint b,
                                       const double delta) const
{
  uint rejectState = mv::fail_state::NO_FAIL;
  double newVolume = volume[b] + delta;
  newDim = *this;
  newDim.SetVolume(b, newVolume);

  //If move would shrink any box axis to be less than 2 * rcut, then
  //automatically reject to prevent errors.
  if ((newDim.halfAx.x[b] < rCut || newDim.halfAx.y[b] < rCut ||
       newDim.halfAx.z[b] < rCut)) {
    std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!"
              << std::endl;
    rejectState = mv::fail_state::VOL_TRANS_WOULD_SHRINK_BOX_BELOW_CUTOFF;
  }
  scale = newDim.axis.Get(b) / axis.Get(b);

  return rejectState;
}

uint BoxDimensionsNonOrth::ExchangeVolume(BoxDimensionsNonOrth & newDim,
    XYZ * scale,
    const double transfer) const
{
  uint state = mv::fail_state::NO_FAIL;
  double vTot = GetTotVolume();
  newDim = *this;

  newDim.SetVolume(0, volume[0] + transfer);
  newDim.SetVolume(1, vTot - newDim.volume[0]);

  //If move would shrink any box axis to be less than 2 * rcut, then
  //automatically reject to prevent errors.
  for (uint b = 0; b < BOX_TOTAL && state == mv::fail_state::NO_FAIL; b++) {
    scale[b] = newDim.axis.Get(b) / axis.Get(b);
    if ((newDim.halfAx.x[b] < rCut || newDim.halfAx.y[b] < rCut ||
         newDim.halfAx.z[b] < rCut)) {
      std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!"
                << std::endl;
      state = state && mv::fail_state::VOL_TRANS_WOULD_SHRINK_BOX_BELOW_CUTOFF;
    }
  }
  return state;
}


void BoxDimensionsNonOrth::SetVolume(const uint b, const double vol)
{
  if(constArea) {
    double ratio = vol / volume[b];
    axis.Scale(b, 1.0, 1.0, ratio);
    halfAx.Scale(b, 1.0, 1.0, ratio);
    cellLength.Scale(b, 1.0, 1.0, ratio);
    //Keep a and b same and change c
    cellBasis[b].Scale(2, ratio);
  } else {
    double ratio = pow(vol / volume[b], (1.0 / 3.0));
    axis.Scale(b, ratio);
    halfAx.Scale(b, ratio);
    cellLength.Scale(b, ratio);
    for(uint i = 0; i < 0; i++) {
      cellBasis[b].Scale(i, ratio);
    }
  }
  volume[b] = vol;
  volInv[b] = 1.0 / volume[b];
  //Calculate new cell dimension
  CalcCellDimensions(b);
}

XYZ BoxDimensionsNonOrth::MinImage(XYZ rawVecRef, const uint b) const
{
  XYZ rawVec = TransformUnSlant(rawVecRef, b);
  rawVecRef = BoxDimensions:: MinImage(rawVec, b);
  rawVecRef = TransformSlant(rawVecRef, b);
  return rawVecRef;
}

void BoxDimensionsNonOrth::WrapPBC(double &x, double &y, double &z,
                                   const uint b) const
{
  //convert XYZ to unslant
  XYZ unwrap(x, y, z);
  XYZ unslant = TransformUnSlant(unwrap, b);
  BoxDimensions::WrapPBC(unslant.x, axis.x[b]);
  BoxDimensions::WrapPBC(unslant.y, axis.y[b]);
  BoxDimensions::WrapPBC(unslant.z, axis.z[b]);
  //convert XYZ to slant
  XYZ slant = TransformSlant(unslant, b);
  x = slant.x;
  y = slant.y;
  z = slant.z;
}

void BoxDimensionsNonOrth::UnwrapPBC(double & x, double & y, double & z,
                                     const uint b, XYZ const& ref) const
{
  //convert XYZ to unslant
  XYZ wrap(x, y, z);
  XYZ unslant = TransformUnSlant(wrap, b);
  XYZ unslantRef = TransformUnSlant(ref, b);
  BoxDimensions::UnwrapPBC(unslant.x, unslantRef.x, axis.x[b], halfAx.x[b]);
  BoxDimensions::UnwrapPBC(unslant.y, unslantRef.y, axis.y[b], halfAx.y[b]);
  BoxDimensions::UnwrapPBC(unslant.z, unslantRef.z, axis.z[b], halfAx.z[b]);
  XYZ unwrap(x, y, z);
  //convert XYZ to slant
  XYZ slant = TransformSlant(unslant, b);
  x = slant.x;
  y = slant.y;
  z = slant.z;
}

//Calculate AxB product
XYZ BoxDimensionsNonOrth::CrossProduct(const XYZ &A, const XYZ &B) const
{
  XYZ temp;
  temp.x = A.y * B.z - A.z * B.y;
  temp.y = A.z * B.x - A.x * B.z;
  temp.z = A.x * B.y - A.y * B.x;

  return temp;
}
