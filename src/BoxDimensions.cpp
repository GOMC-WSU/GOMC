/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
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


    axis.Set(b, cellBasis[b].Length(0), cellBasis[b].Length(1),
             cellBasis[b].Length(2));

    if(axis.Get(b).Min() < 2.0 * rCut[b]) {
      printf("Error: Cutoff value is larger than half of minimum BOX%d length!\n", b);
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
                                const uint b, const double delta) const
{
  uint rejectState = mv::fail_state::NO_FAIL;
  double newVolume = volume[b] + delta;
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

double BoxDimensions::GetTotVolume(const uint b1, const uint b2) const
{
  return (volume[b1] + volume[b2]);
}

void BoxDimensions::SetVolume(const uint b, const double vol)
{
  if(constArea) {
    double ratio = vol / volume[b];
    axis.Scale(b, 1.0, 1.0, ratio);
    halfAx.Scale(b, 1.0, 1.0, ratio);
  } else {
    double ratio = cbrt(vol / volume[b]);
    axis.Scale(b, ratio);
    halfAx.Scale(b, ratio);
  }
  volume[b] = vol;
  volInv[b] = 1.0 / vol;
}


void BoxDimensions::WrapPBC(double &x, double &y, double &z, const uint b) const
{
  WrapPBC(x, axis.x[b]);
  WrapPBC(y, axis.y[b]);
  WrapPBC(z, axis.z[b]);
}

void BoxDimensions::WrapPBC(double &x, double &y, double &z, const uint b,
                      const bool &pbcX, const bool &pbcY, const bool &pbcZ) const
{
  if(pbcX) {WrapPBC(x, axis.x[b]);}
  if(pbcY) {WrapPBC(y, axis.y[b]);}
  if(pbcZ) {WrapPBC(z, axis.z[b]);}
}

void BoxDimensions::UnwrapPBC(double & x, double & y, double & z,
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
  if ( v >= ax ) //if +, wrap out to low end, on boundary will wrap to zero
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
  rawVec.x = MinImageSigned(rawVec.x, axis.x[b], halfAx.x[b]);
  rawVec.y = MinImageSigned(rawVec.y, axis.y[b], halfAx.y[b]);
  rawVec.z = MinImageSigned(rawVec.z, axis.z[b], halfAx.z[b]);
  return rawVec;
}

XYZ BoxDimensions::MinImage_X(XYZ rawVec, const uint b) const
{
  rawVec.x = MinImageSigned(rawVec.x, axis.x[b], halfAx.x[b]);
  return rawVec;
}

XYZ BoxDimensions::MinImage_Y(XYZ rawVec, const uint b) const
{
  rawVec.y = MinImageSigned(rawVec.y, axis.y[b], halfAx.y[b]);
  return rawVec;
}

XYZ BoxDimensions::MinImage_Z(XYZ rawVec, const uint b) const
{
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
                                                                      + _________
 NUMBER_OF_ATTRIBUTES                                                 // 26 / 51
  */

#if GOMC_LIB_MPI
  std::vector<double> BoxDimensions::SerializeBoxDimObject(){
    #if ENSEMBLE == GEMC
      int NUMBER_OF_ATTRIBUTES = 51;
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
    #else
      int NUMBER_OF_ATTRIBUTES = 26;
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
    #endif

    return serialBoxDim;
  }

    void BoxDimensions::ReadFromSerializedBoxDimObject(std::vector<double> & serialBoxDim){
    
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

    #endif

  }
#endif