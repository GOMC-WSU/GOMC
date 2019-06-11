/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef COORDINATES_H
#define COORDINATES_H

#include "BasicTypes.h"
#include "Molecules.h" //For start
#include "BoxDimensions.h" //For pbc wrapping
#include "BoxDimensionsNonOrth.h"
#include "XYZArray.h" //Parent class
#include "MoleculeLookup.h" //For box iterators used in initial assignment
#include "COM.h"
#include "PRNG.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>

//Coordinates array
class Coordinates : public XYZArray
{
public:
  //Declare a set of coordinates with no data (but must have proper frame
  //of reference).
  Coordinates(BoxDimensions & box, COM & com,
              MoleculeLookup & molLook, PRNG & prng, Molecules const& mol) :
    boxDimRef(box), comRef(com), prngRef(prng), molLookRef(molLook),
    molRef(mol) {}

  Coordinates& operator=(Coordinates const& rhs)
  {
    this->XYZArray::operator=(rhs);
    return *this;
  }

  //Init from the coordinates grabbed from pdb file read.
  void InitFromPDB(pdb_setup::Atoms const& atoms);

  // to see if they are within defined volume or not.
  void CheckCoordinate();

  //Translate by a random amount
  void TranslateRand(XYZArray & dest, XYZ & newCOM, uint & pStart,
                     uint & pLen, const uint m, const uint b,
                     real max);

  //Rotate by a random amount.
  void RotateRand(XYZArray & dest,  uint & pStart, uint & pLen, const uint m,
                  const uint b, const real max);

  //scale all in each mol newCOM[m]/oldCOM[m]
  void VolumeTransferTranslate
  (uint & state, Coordinates &dest, COM & newCOM, BoxDimensions & newDim,
   COM const& oldCOM, const real max, const uint *box) const;

  //Helper for TranslateAll
  void TranslateOneBox(Coordinates & dest, COM & newCOM, COM const& oldCOM,
                       BoxDimensions const& newDim, const uint b,
                       const XYZ& scale) const;

private:

  BoxDimensions & boxDimRef;
  COM & comRef;
  PRNG & prngRef;
  MoleculeLookup & molLookRef;
  Molecules const& molRef;
};



#endif /*COORDINATES_H*/
