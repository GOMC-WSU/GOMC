/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef BOX_DIMENSIONS_NONORTHO_H
#define BOX_DIMENSIONS_NONORTHO_H

#include "BoxDimensions.h"

class BoxDimensionsNonOrth : public BoxDimensions {
public:
  BoxDimensionsNonOrth() : BoxDimensions() {
    cellLength.Init(BOX_TOTAL);
    for (uint b = 0; b < BOX_TOTAL; b++) {
      cellBasis[b] = XYZArray(3);
      cellBasis_Inv[b] = XYZArray(3);
    }
  }
  BoxDimensionsNonOrth(BoxDimensionsNonOrth const &other)
      : BoxDimensions(other) {
    cellLength.Init(BOX_TOTAL);
    other.cellLength.CopyRange(cellLength, 0, 0, BOX_TOTAL);
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      cellBasis_Inv[b] = XYZArray(3);
      other.cellBasis_Inv[b].CopyRange(cellBasis_Inv[b], 0, 0, 3);
    }
  }

  ~BoxDimensionsNonOrth(){};

  BoxDimensionsNonOrth &operator=(BoxDimensionsNonOrth const &other);
  bool operator==(BoxDimensionsNonOrth const &other);

  void Init(config_setup::RestartSettings const &restart,
            config_setup::Volume const &confVolume,
            pdb_setup::Cryst1 const &cryst, Forcefield const &ff) override;

  void SetVolume(const uint b, const double vol) override;

  uint ShiftVolume(BoxDimensionsNonOrth &newDim, XYZ &scale, const uint b,
                   const double delta) const;

  //! Calculate and execute volume exchange based on transfer
  uint ExchangeVolume(BoxDimensionsNonOrth &newDim, XYZ *scale,
                      const double transfer, const uint *box) const;

  // Construct cell basis based on new axis dimension
  void CalcCellDimensions(const uint b);

  // Vector btwn two points, accounting for PBC, on an individual axis
  XYZ MinImage(XYZ rawVecRef, const uint b) const override;

  // Apply PBC, on X axis
  XYZ MinImage_X(XYZ rawVec, const uint b) const override;
  // Apply PBC, on Y axis
  XYZ MinImage_Y(XYZ rawVec, const uint b) const override;
  // Apply PBC, on Z axis
  XYZ MinImage_Z(XYZ rawVec, const uint b) const override;

  // Wrap one coordinate
  void WrapPBC(double &x, double &y, double &z, const uint b) const override;

  // Wrap one coordinate and check for PBC
  void WrapPBC(double &x, double &y, double &z, const uint b, const bool &pbcX,
               const bool &pbcY, const bool &pbcZ) const override;

  // Unwrap one coordinate
  void UnwrapPBC(double &x, double &y, double &z, const uint b,
                 XYZ const &ref) const override;

  // Transform A to unslant coordinate
  XYZ TransformUnSlant(const XYZ &A, const uint b) const override;

  // Transform A to slant coordinate
  XYZ TransformSlant(const XYZ &A, const uint b) const override;

  // private:
  XYZArray cellBasis_Inv[BOX_TOTAL]; // inverse cell matrix for each box
  XYZArray cellLength;               // Length of a, b, c for each box
};

// Calculate inverse transform
inline XYZ BoxDimensionsNonOrth::TransformUnSlant(const XYZ &A,
                                                  const uint b) const {
  XYZ temp;

  temp.x = A.x * cellBasis_Inv[b].Get(0).x + A.y * cellBasis_Inv[b].Get(1).x +
           A.z * cellBasis_Inv[b].Get(2).x;
  temp.y = A.x * cellBasis_Inv[b].Get(0).y + A.y * cellBasis_Inv[b].Get(1).y +
           A.z * cellBasis_Inv[b].Get(2).y;
  temp.z = A.x * cellBasis_Inv[b].Get(0).z + A.y * cellBasis_Inv[b].Get(1).z +
           A.z * cellBasis_Inv[b].Get(2).z;
  return temp;
}

// Calculate transform
inline XYZ BoxDimensionsNonOrth::TransformSlant(const XYZ &A,
                                                const uint b) const {
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
