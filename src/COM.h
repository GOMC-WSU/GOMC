/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef COM_H
#define COM_H

#include "BasicTypes.h"
#include "Molecules.h" //For start
#include "BoxDimensions.h" //For pbc wrapping
#include "BoxDimensionsNonOrth.h"
#include "XYZArray.h" //Parent class
#include "MoleculeLookup.h" //For box iterators used in initial assignment

//COM array
class COM : public XYZArray
{
public:

  //void operator=(const COM& b)
  COM& operator=(const COM& b)
  {
    this->XYZArray::operator=(b);
    return *this;
  }
  //Declare a set of coordinates with no data (but must have proper frame
  //of reference).
  COM(BoxDimensions & box, XYZArray & coordinates, MoleculeLookup & molLook,
      Molecules const& mol) :
    boxDimRef(box), coordRef(coordinates), molLookRef(molLook), molRef(mol)
  {
  }

  //Init from the coordinates grabbed from pdb file read.
  void CalcCOM();
  void SetNew(const uint m, const uint b);

private:

  BoxDimensions & boxDimRef;
  XYZArray & coordRef;
  MoleculeLookup & molLookRef;
  Molecules const& molRef;
};

inline void COM::CalcCOM()
{
  MoleculeLookup::box_iterator current, end;
  uint pStart = 0, pStop = 0, pLen = 0;
  XYZArray::Init(molRef.count);
  for (uint b = 0; b < BOX_TOTAL; b++) {
    current = molLookRef.BoxBegin(b);
    end = molLookRef.BoxEnd(b);
    while (current != end) {
      molRef.GetRange(pStart, pStop, pLen, *current);
      boxDimRef.UnwrapPBC(coordRef, pStart, pStop,
                          b, coordRef.Get(pStart));
      Set(*current, 0, 0, 0);
      for (uint p = pStart; p < pStop; p++)
        Add(coordRef, *current, p);
      boxDimRef.WrapPBC(coordRef, pStart, pStop, b);
      Scale(*current, 1.0 / (real)(pLen));
      boxDimRef.WrapPBC(x[*current], y[*current], z[*current], b);
      ++current;
    }
  }
}

inline void COM::SetNew(const uint m, const uint b)
{
  uint pStart = 0, pStop = 0, pLen = 0;
  molRef.GetRange(pStart, pStop, pLen, m);
  boxDimRef.UnwrapPBC(coordRef, pStart, pStop,
                      b, coordRef.Get(pStart));
  Set(m, 0.0, 0.0, 0.0);
  for (uint p = pStart; p < pStop; p++)
    Add(coordRef, m, p);
  boxDimRef.WrapPBC(coordRef, pStart, pStop, b);
  Scale(m, 1.0 / (real)(pLen));
  boxDimRef.WrapPBC(x[m], y[m], z[m], b);
}

#endif /*COM_H*/
