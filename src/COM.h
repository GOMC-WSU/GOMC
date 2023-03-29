/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef COM_H
#define COM_H

#include "BasicTypes.h"
#include "BoxDimensions.h" //For pbc wrapping
#include "BoxDimensionsNonOrth.h"
#include "MoleculeLookup.h" //For box iterators used in initial assignment
#include "Molecules.h"      //For start
#include "XYZArray.h"       //Parent class

// COM array
class COM : public XYZArray {
public:
  // void operator=(const COM& b)
  COM &operator=(const COM &b) {
    this->XYZArray::operator=(b);
    return *this;
  }
  // Declare a set of coordinates with no data (but must have proper frame
  // of reference).
  COM(BoxDimensions &box, XYZArray &coordinates, MoleculeLookup &molLook,
      Molecules const &mol)
      : boxDimRef(box), coordRef(coordinates), molLookRef(molLook),
        molRef(mol) {}

  // Init from the coordinates grabbed from pdb file read.
  void CalcCOM();
  void SetNew(const uint m, const uint b);

private:
  BoxDimensions &boxDimRef;
  XYZArray &coordRef;
  MoleculeLookup &molLookRef;
  Molecules const &molRef;
};

inline void COM::CalcCOM() {
  MoleculeLookup::box_iterator current, end;
  uint pStart = 0, pStop = 0, pLen = 0;
  XYZArray::Init(molRef.count);
  for (uint b = 0; b < BOX_TOTAL; b++) {
    current = molLookRef.BoxBegin(b);
    end = molLookRef.BoxEnd(b);
    while (current != end) {
      molRef.GetRange(pStart, pStop, pLen, *current);
      boxDimRef.UnwrapPBC(coordRef, pStart, pStop, b, coordRef.Get(pStart));
      Set(*current, 0, 0, 0);
      for (uint p = pStart; p < pStop; p++)
        Add(coordRef, *current, p);
      boxDimRef.WrapPBC(coordRef, pStart, pStop, b);
      Scale(*current, 1.0 / (double)(pLen));
      boxDimRef.WrapPBC(x[*current], y[*current], z[*current], b);
      ++current;
    }
  }
}

inline void COM::SetNew(const uint moleculeIndex, const uint box) {
  uint atomStart = 0, atomStop = 0, atomLen = 0;

  // get the range of atoms for molecueIndex
  molRef.GetRange(atomStart, atomStop, atomLen, moleculeIndex);

  // unwrap molecule
  boxDimRef.UnwrapPBC(coordRef, atomStart, atomStop, box,
                      coordRef.Get(atomStart));

  // set COM of moleculeIndex to all 0.0
  Set(moleculeIndex, 0.0, 0.0, 0.0);

  // add the coordinates of all atoms to COM of molecule moleculeIndex
  for (uint currentAtom = atomStart; currentAtom < atomStop; currentAtom++)
    Add(coordRef, moleculeIndex, currentAtom);

  // wrap atoms again
  boxDimRef.WrapPBC(coordRef, atomStart, atomStop, box);

  // to take average (calculate COM) divide the sum by atom length
  Scale(moleculeIndex, 1.0 / (double)(atomLen));

  // wrap the COM
  boxDimRef.WrapPBC(x[moleculeIndex], y[moleculeIndex], z[moleculeIndex], box);
}

#endif /*COM_H*/
