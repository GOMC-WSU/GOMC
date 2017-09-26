/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.1
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BOX_DIMENSIONS_NONORTHO_H
#define BOX_DIMENSIONS_NONORTHO_H

#include "BoxDimensions.h"


class BoxDimensionsNonOrth : public BoxDimensions
{
public:
  BoxDimensionsNonOrth() {}
  BoxDimensionsNonOrth(BoxDimensionsNonOrth const& other) : BoxDimensions(other)
  {
    other.cellLength.CopyRange(cellLength, 0, 0, BOX_TOTAL);
    other.faceLength.CopyRange(faceLength, 0, 0, BOX_TOTAL);
    for (uint b = 0; b < BOX_TOTAL; ++b)
    {   
      other.cellBasis_Inv[b].CopyRange(cellBasis_Inv[b], 0, 0, 3);
    }
  }

  virtual BoxDimensionsNonOrth& operator=(BoxDimensionsNonOrth const& other);

  virtual void Init(config_setup::RestartSettings const& restart,
             config_setup::Volume const& confVolume,
             pdb_setup::Cryst1 const& cryst, double rc, double rcSq);

  virtual void SetVolume(const uint b, const double vol);

  //Construct cell basis based on new axis dimension
  void CalcCellDimensions();

  //Vector btwn two points, accounting for PBC, on an individual axis
  virtual XYZ MinImage(XYZ rawVec, const uint b) const;

  //Unwrap one coordinate.
  virtual void WrapPBC(double &x, double &y, double &z, const uint b) const;

  //Unwrap one coordinate.
  virtual void UnwrapPBC(double & x, double & y, double & z,
			 const uint b, XYZ const& ref) const;

  //Transform A to unslant coordinate
  XYZ TransformUnSlant(const XYZ &A, const uint b) const;

  //Transform A to slant coordinate
  XYZ TransformSlant(const XYZ &A, const uint b) const;

//private:
  XYZArray cellBasis_Inv[BOX_TOTAL]; //inverse cell matrix for each box
  XYZArray cellLength;                //Length of a, b, c for each box
  XYZArray faceLength;                //Length between two faces for each box

  XYZ CrossProduct(const XYZ &A, const XYZ &B) const;   //Calc AxB product  
};



#endif /*BOX_DIMENSIONS_NONORTHO_H*/
