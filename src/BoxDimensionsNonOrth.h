/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BOX_DIMENSIONS_NONORTHO_H
#define BOX_DIMENSIONS_NONORTHO_H

#include "BoxDimensions.h"


class BoxDimensionsNonOrth : public BoxDimensions
{
public:
  BoxDimensionsNonOrth() : BoxDimensions()
  {
    cellLength.Init(BOX_TOTAL);
    for (uint b = 0; b < BOX_TOTAL; b++) {
      cellBasis[b] = XYZArray(3);
      cellBasis_Inv[b] = XYZArray(3);
    }
  }
  BoxDimensionsNonOrth(BoxDimensionsNonOrth const& other) : BoxDimensions(other)
  {
    cellLength.Init(BOX_TOTAL);
    other.cellLength.CopyRange(cellLength, 0, 0, BOX_TOTAL);
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      cellBasis_Inv[b] = XYZArray(3);
      other.cellBasis_Inv[b].CopyRange(cellBasis_Inv[b], 0, 0, 3);
    }
  }

  virtual BoxDimensionsNonOrth& operator=(BoxDimensionsNonOrth const& other);

  virtual void Init(config_setup::RestartSettings const& restart,
                    config_setup::Volume const& confVolume,
                    pdb_setup::Cryst1 const& cryst,
                    Forcefield const &ff);

  virtual void SetVolume(const uint b, const real vol);

  virtual uint ShiftVolume(BoxDimensionsNonOrth & newDim, XYZ & scale,
                           const uint b, const real delta) const;

  //!Calculate and execute volume exchange based on transfer
  virtual uint ExchangeVolume(BoxDimensionsNonOrth & newDim, XYZ * scale,
                              const real transfer, const uint *box) const;

  //Construct cell basis based on new axis dimension
  void CalcCellDimensions(const uint b);

  //Vector btwn two points, accounting for PBC, on an individual axis
  virtual XYZ MinImage(XYZ rawVecRef, const uint b) const;

  //Unwrap one coordinate.
  virtual void WrapPBC(real &x, real &y, real &z, const uint b) const;

  //Unwrap one coordinate.
  virtual void UnwrapPBC(real & x, real & y, real & z,
                         const uint b, XYZ const& ref) const;

  //Transform A to unslant coordinate
  XYZ TransformUnSlant(const XYZ &A, const uint b) const;

  //Transform A to slant coordinate
  XYZ TransformSlant(const XYZ &A, const uint b) const;

//private:
  XYZArray cellBasis_Inv[BOX_TOTAL]; //inverse cell matrix for each box
  XYZArray cellLength;                //Length of a, b, c for each box
};

//Calculate transform
inline XYZ BoxDimensionsNonOrth::TransformUnSlant(const XYZ &A,
    const uint b) const
{
  XYZ temp;

  temp.x = A.x * cellBasis_Inv[b].Get(0).x + A.y * cellBasis_Inv[b].Get(1).x +
           A.z * cellBasis_Inv[b].Get(2).x;
  temp.y = A.x * cellBasis_Inv[b].Get(0).y + A.y * cellBasis_Inv[b].Get(1).y +
           A.z * cellBasis_Inv[b].Get(2).y;
  temp.z = A.x * cellBasis_Inv[b].Get(0).z + A.y * cellBasis_Inv[b].Get(1).z +
           A.z * cellBasis_Inv[b].Get(2).z;
  return temp;
}

//Calculate transform
inline XYZ BoxDimensionsNonOrth::TransformSlant(const XYZ &A, const uint b) const
{
  XYZ temp;

  temp.x = A.x * cellBasis[b].Get(0).x + A.y * cellBasis[b].Get(1).x +
           A.z * cellBasis[b].Get(2).x;
  temp.y = A.x * cellBasis[b].Get(0).y + A.y * cellBasis[b].Get(1).y +
           A.z * cellBasis[b].Get(2).y;
  temp.z = A.x * cellBasis[b].Get(0).z + A.y * cellBasis[b].Get(1).z +
           A.z * cellBasis[b].Get(2).z;
  return temp;
}

#endif /*BOX_DIMENSIONS_NONORTHO_H*/
