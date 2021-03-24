/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "BoxDimensionsNonOrth.h"
#include "BoxDimensions.h"
#include "GeomLib.h"
#include "MoveConst.h" //For cutoff-related fail condition

using namespace geom;

void BoxDimensionsNonOrth::Init(config_setup::RestartSettings const& restart,
                                config_setup::Volume const& confVolume,
                                pdb_setup::Cryst1 const& cryst,
                                Forcefield const &ff)
{
  for (uint b = 0; b < BOX_TOTAL; b++) {
    rCut[b] = std::max(ff.rCut, ff.rCutCoulomb[b]);
    rCutSq[b] = rCut[b] * rCut[b];
    minVol[b] = 8.0 * rCutSq[b] * rCut[b] + 0.001;
    if(restart.enable && cryst.hasVolume[b]) {
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
      double cosBSq = beta * beta;
      double cosGSq = gamma * gamma;
      double temp = (alpha - beta * gamma) / (sqrt(1.0 - cosGSq));
      cellBasis[b].Set(0, 1.0, 0.0, 0.0);
      cellBasis[b].Set(1, gamma, sqrt(1.0 - cosGSq), 0.0);
      cellBasis[b].Set(2, beta, temp, sqrt(1.0 - cosBSq - temp * temp));
      cellBasis[b].Scale(0, axis.Get(b).x);
      cellBasis[b].Scale(1, axis.Get(b).y);
      cellBasis[b].Scale(2, axis.Get(b).z);
    } else if (restart.enable && cryst.hasCellBasis[b]) {
      cryst.cellBasis[b].CopyRange(cellBasis[b], 0, 0, 3);
    } else if(confVolume.hasVolume) {
      confVolume.axis[b].CopyRange(cellBasis[b], 0, 0, 3);
    } else {
      fprintf(stderr,
              "Error: Cell Basis not specified in XSC, PDB, or in.conf files.\n");
      exit(EXIT_FAILURE);
    }

    //Print Box dimension info
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

    //Find the length of a, b, c
    cellLength.Set(b, cellBasis[b].Length(0), cellBasis[b].Length(1),
                   cellBasis[b].Length(2));
    //Find Cosine Angle of alpha, beta and gamma
    cosAngle[b][0] = Dot(cellBasis[b].Get(1), cellBasis[b].Get(2)) /
                     (cellLength.Get(b).y * cellLength.Get(b).z);
    cosAngle[b][1] = Dot(cellBasis[b].Get(0), cellBasis[b].Get(2)) /
                     (cellLength.Get(b).x * cellLength.Get(b).z);
    cosAngle[b][2] = Dot(cellBasis[b].Get(0), cellBasis[b].Get(1)) /
                     (cellLength.Get(b).x * cellLength.Get(b).y);
    //Calculate Cross Product
    XYZ bxc = Cross(cellBasis[b].Get(1), cellBasis[b].Get(2));
    //Calculate volume = A.(B x C)
    volume[b] = std::abs(Dot(cellBasis[b].Get(0), bxc));
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

    if(axis.Get(b).Min() < 2.0 * rCut[b]) {
      printf("Error: Cutoff value is larger than half of minimum BOX%d length!\n", b);
      exit(EXIT_FAILURE);
    }
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
  other.cellLength.CopyRange(cellLength, 0, 0, BOX_TOTAL);
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

uint BoxDimensionsNonOrth::ExchangeVolume(BoxDimensionsNonOrth & newDim,
    XYZ * scale, const double transfer, const uint *box) const
{
  uint state = mv::fail_state::NO_FAIL;
  double vTot = GetTotVolume(box[0], box[1]);
  newDim = *this;

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


void BoxDimensionsNonOrth::SetVolume(const uint b, const double vol)
{
  if(constArea) {
    double ratio = vol / volume[b];
    axis.Scale(b, 1.0, 1.0, ratio);
    halfAx.Scale(b, 1.0, 1.0, ratio);
    cellLength.Scale(b, 1.0, 1.0, ratio);
  } else {
    double ratio = cbrt(vol / volume[b]);
    axis.Scale(b, ratio);
    halfAx.Scale(b, ratio);
    cellLength.Scale(b, ratio);
  }
  volume[b] = vol;
  volInv[b] = 1.0 / volume[b];
  //Calculate new cell dimension
  CalcCellDimensions(b);
}

XYZ BoxDimensionsNonOrth::MinImage(XYZ rawVecRef, const uint b) const
{
  XYZ rawVec = TransformUnSlant(rawVecRef, b);
  rawVecRef = BoxDimensions::MinImage(rawVec, b);
  rawVecRef = TransformSlant(rawVecRef, b);
  return rawVecRef;
}

XYZ BoxDimensionsNonOrth::MinImage_X(XYZ rawVecRef, const uint b) const
{
  XYZ rawVec = TransformUnSlant(rawVecRef, b);
  rawVecRef = BoxDimensions::MinImage_X(rawVec, b);
  rawVecRef = TransformSlant(rawVecRef, b);
  return rawVecRef;
}

XYZ BoxDimensionsNonOrth::MinImage_Y(XYZ rawVecRef, const uint b) const
{
  XYZ rawVec = TransformUnSlant(rawVecRef, b);
  rawVecRef = BoxDimensions::MinImage_Y(rawVec, b);
  rawVecRef = TransformSlant(rawVecRef, b);
  return rawVecRef;
}

XYZ BoxDimensionsNonOrth::MinImage_Z(XYZ rawVecRef, const uint b) const
{
  XYZ rawVec = TransformUnSlant(rawVecRef, b);
  rawVecRef = BoxDimensions::MinImage_Z(rawVec, b);
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

void BoxDimensionsNonOrth::WrapPBC(double &x, double &y, double &z,
                                   const uint b, const bool &pbcX,
                                   const bool &pbcY, const bool &pbcZ) const
{
  //convert XYZ to unslant
  XYZ unwrap(x, y, z);
  XYZ unslant = TransformUnSlant(unwrap, b);
  if(pbcX) {BoxDimensions::WrapPBC(unslant.x, axis.x[b]);}
  if(pbcY) {BoxDimensions::WrapPBC(unslant.y, axis.y[b]);}
  if(pbcZ) {BoxDimensions::WrapPBC(unslant.z, axis.z[b]);}
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
  //convert XYZ to slant
  XYZ slant = TransformSlant(unslant, b);
  x = slant.x;
  y = slant.y;
  z = slant.z;
}

/*
                                                                  #BOXES  1 / 2
                                                                        _________
  XYZArray axis;                  //x, y, z dimensions of each box (a)    3 / 6
  XYZArray halfAx;               //x, y, z dimensions / 2 of each box (a) 3 / 6
  XYZArray cellBasis[BOX_TOTAL];  //x, y, z vector, 3 for each box        9 / 18
  double volume[BOX_TOTAL];       //volume of each box in (a^3)           1 / 2
  double volInv[BOX_TOTAL];       //inverse volume of each box in (a^-3)  1 / 2
  double cosAngle[BOX_TOTAL][3];  //alpha, beta, gamma for each box       3 / 6
  double rCut[BOX_TOTAL];                                               //1 / 2
  double rCutSq[BOX_TOTAL];                                             //1 / 2
  double minVol[BOX_TOTAL];                                             //1 / 2
  bool cubic[BOX_TOTAL], orthogonal[BOX_TOTAL], constArea;              //3 / 5

// NONORTH
  XYZArray cellBasis_Inv[BOX_TOTAL]; //inverse cell matrix for each box   9 / 18
  XYZArray cellLength;                //Length of a, b, c for each box    3 / 6

                                                                      + _________
 NUMBER_OF_ATTRIBUTES                                                 // 38 / 75
  */

#if GOMC_LIB_MPI
  std::vector<double> BoxDimensionsNonOrth::SerializeBoxDimObject(){
    #if ENSEMBLE == GEMC
      int NUMBER_OF_ATTRIBUTES = 75;
      std::vector<double> serialBoxDim(NUMBER_OF_ATTRIBUTES);
      serialBoxDim[0] = axis[0].x;
      serialBoxDim[1] = axis[0].y;
      serialBoxDim[2] = axis[0].z;
      serialBoxDim[3] = axis[1].x;
      serialBoxDim[4] = axis[1].y;
      serialBoxDim[5] = axis[1].z;
      serialBoxDim[6] = halfAx[0].x;
      serialBoxDim[7] = halfAx[0].y;
      serialBoxDim[8] = halfAx[0].z;
      serialBoxDim[9] = halfAx[1].x;
      serialBoxDim[10] = halfAx[1].y;
      serialBoxDim[11] = halfAx[1].z;
      memcpy(&serialBoxDim[12], cellBasis[0].x, 3*sizeof(double));
      memcpy(&serialBoxDim[15], cellBasis[0].y, 3*sizeof(double));
      memcpy(&serialBoxDim[18], cellBasis[0].z, 3*sizeof(double));
      memcpy(&serialBoxDim[21], cellBasis[1].x, 3*sizeof(double));
      memcpy(&serialBoxDim[24], cellBasis[1].y, 3*sizeof(double));
      memcpy(&serialBoxDim[27], cellBasis[1].z, 3*sizeof(double));
      serialBoxDim[30] = volume[0];
      serialBoxDim[31] = volume[1];
      serialBoxDim[32] = volInv[0];
      serialBoxDim[33] = volInv[1];
      memcpy(&serialBoxDim[34], cosAngle[0], 3*sizeof(double));
      memcpy(&serialBoxDim[37], cosAngle[1], 3*sizeof(double));
      serialBoxDim[40] = rCut[0];
      serialBoxDim[41] = rCut[1];
      serialBoxDim[42] = rCutSq[0];
      serialBoxDim[43] = rCutSq[1];
      serialBoxDim[44] = minVol[0];
      serialBoxDim[45] = minVol[1];
      serialBoxDim[46] = (double)cubic[0];
      serialBoxDim[47] = (double)cubic[1];
      serialBoxDim[48] = (double)orthogonal[0];
      serialBoxDim[49] = (double)orthogonal[1];
      serialBoxDim[50] = (double)constArea;

      memcpy(&serialBoxDim[51], cellBasis_Inv[0].x, 3*sizeof(double));
      memcpy(&serialBoxDim[54], cellBasis_Inv[0].y, 3*sizeof(double));
      memcpy(&serialBoxDim[57], cellBasis_Inv[0].z, 3*sizeof(double));
      memcpy(&serialBoxDim[60], cellBasis_Inv[1].x, 3*sizeof(double));
      memcpy(&serialBoxDim[63], cellBasis_Inv[1].y, 3*sizeof(double));
      memcpy(&serialBoxDim[66], cellBasis_Inv[1].z, 3*sizeof(double));

      serialBoxDim[69] = cellLength[0].x;
      serialBoxDim[70] = cellLength[0].y;
      serialBoxDim[71] = cellLength[0].z;
      serialBoxDim[72] = cellLength[1].x;
      serialBoxDim[73] = cellLength[1].y;
      serialBoxDim[74] = cellLength[1].z;

    #else
      int NUMBER_OF_ATTRIBUTES = 38;
      std::vector<double> serialBoxDim(NUMBER_OF_ATTRIBUTES);
      serialBoxDim[0] = axis[0].x;
      serialBoxDim[1] = axis[0].y;
      serialBoxDim[2] = axis[0].z;
      serialBoxDim[3] = halfAx[0].x;
      serialBoxDim[4] = halfAx[0].y;
      serialBoxDim[5] = halfAx[0].z;
      memcpy(&serialBoxDim[6], cellBasis[0].x, 3*sizeof(double));
      memcpy(&serialBoxDim[9], cellBasis[0].y, 3*sizeof(double));
      memcpy(&serialBoxDim[12], cellBasis[0].z, 3*sizeof(double));
      serialBoxDim[15] = volume[0];
      serialBoxDim[16] = volInv[0];
      memcpy(&serialBoxDim[17], cosAngle[0], 3*sizeof(double));
      serialBoxDim[20] = rCut[0];
      serialBoxDim[21] = rCutSq[0];
      serialBoxDim[22] = minVol[0];
      serialBoxDim[23] = (double)cubic[0];
      serialBoxDim[24] = (double)orthogonal[0];
      serialBoxDim[25] = (double)constArea;

      memcpy(&serialBoxDim[26], cellBasis_Inv[0].x, 3*sizeof(double));
      memcpy(&serialBoxDim[29], cellBasis_Inv[0].y, 3*sizeof(double));
      memcpy(&serialBoxDim[32], cellBasis_Inv[0].z, 3*sizeof(double));

      serialBoxDim[35] = cellLength[0].x;
      serialBoxDim[36] = cellLength[0].y;
      serialBoxDim[37] = cellLength[0].z;

    #endif

    return serialBoxDim;
  }

    void BoxDimensionsNonOrth::ReadFromSerializedBoxDimObject(std::vector<double> & serialBoxDim){
    
    #if ENSEMBLE == GEMC

      axis[0].x = serialBoxDim[0];
      axis[0].y = serialBoxDim[1];
      axis[0].z = serialBoxDim[2];
      axis[1].x = serialBoxDim[3];
      axis[1].y = serialBoxDim[4];
      axis[1].z = serialBoxDim[5];
      halfAx[0].x = serialBoxDim[6];
      halfAx[0].y = serialBoxDim[7];
      halfAx[0].z = serialBoxDim[8];
      halfAx[1].x = serialBoxDim[9];
      halfAx[1].y = serialBoxDim[10];
      halfAx[1].z = serialBoxDim[11];

      memcpy(cellBasis[0].x, &serialBoxDim[12], 3*sizeof(double));
      memcpy(cellBasis[0].y, &serialBoxDim[15], 3*sizeof(double));
      memcpy(cellBasis[0].z, &serialBoxDim[18], 3*sizeof(double));
      memcpy(cellBasis[1].x, &serialBoxDim[21], 3*sizeof(double));
      memcpy(cellBasis[1].y, &serialBoxDim[24], 3*sizeof(double));
      memcpy(cellBasis[1].z, &serialBoxDim[27], 3*sizeof(double));

      volume[0] = serialBoxDim[30];
      volume[1] = serialBoxDim[31];
      volInv[0] = serialBoxDim[32];
      volInv[1] = serialBoxDim[33];

      memcpy(cosAngle[0], &serialBoxDim[34], 3*sizeof(double));
      memcpy(cosAngle[1], &serialBoxDim[37], 3*sizeof(double));

      rCut[0] = serialBoxDim[40];
      rCut[1] = serialBoxDim[41];
      rCutSq[0] = serialBoxDim[42];
      rCutSq[1] = serialBoxDim[43];
      minVol[0] = serialBoxDim[44];
      minVol[1] = serialBoxDim[45];
      cubic[0] = (bool)serialBoxDim[46];
      cubic[1] = (bool)serialBoxDim[47];
      orthogonal[0] = (bool)serialBoxDim[48];
      orthogonal[1] = (bool)serialBoxDim[49];
      constArea = (bool)serialBoxDim[50];

      memcpy(cellBasis_Inv[0].x, &serialBoxDim[51], 3*sizeof(double));
      memcpy(cellBasis_Inv[0].y, &serialBoxDim[54], 3*sizeof(double));
      memcpy(cellBasis_Inv[0].z, &serialBoxDim[57], 3*sizeof(double));
      memcpy(cellBasis_Inv[1].x, &serialBoxDim[60], 3*sizeof(double));
      memcpy(cellBasis_Inv[1].y, &serialBoxDim[63], 3*sizeof(double));
      memcpy(cellBasis_Inv[1].z, &serialBoxDim[66], 3*sizeof(double));

      cellLength[0].x = serialBoxDim[69];
      cellLength[0].y = serialBoxDim[70];
      cellLength[0].z = serialBoxDim[71];
      cellLength[1].x = serialBoxDim[72];
      cellLength[1].y = serialBoxDim[73];
      cellLength[1].z = serialBoxDim[74];

    #else

      axis[0].x = serialBoxDim[0];
      axis[0].y = serialBoxDim[1];
      axis[0].z = serialBoxDim[2];
      halfAx[0].x = serialBoxDim[3];
      halfAx[0].y = serialBoxDim[4];
      halfAx[0].z = serialBoxDim[5];

      memcpy(cellBasis[0].x, &serialBoxDim[6], 3*sizeof(double));
      memcpy(cellBasis[0].y, &serialBoxDim[9], 3*sizeof(double));
      memcpy(cellBasis[0].z, &serialBoxDim[12], 3*sizeof(double));

      volume[0] = serialBoxDim[15];
      volInv[0] = serialBoxDim[16];

      memcpy(cosAngle[0], &serialBoxDim[17], 3*sizeof(double));

      rCut[0] = serialBoxDim[20];
      rCutSq[0] = serialBoxDim[21];
      minVol[0] = serialBoxDim[22];
      cubic[0] = (bool)serialBoxDim[23];
      orthogonal[0] = (bool)serialBoxDim[24];
      constArea = (bool)serialBoxDim[25];

      memcpy(cellBasis_Inv[0].x, &serialBoxDim[26], 3*sizeof(double));
      memcpy(cellBasis_Inv[0].y, &serialBoxDim[29], 3*sizeof(double));
      memcpy(cellBasis_Inv[0].z, &serialBoxDim[32], 3*sizeof(double));

      cellLength[0].x = serialBoxDim[35];
      cellLength[0].y = serialBoxDim[36];
      cellLength[0].z = serialBoxDim[37];

    #endif

  }
#endif
