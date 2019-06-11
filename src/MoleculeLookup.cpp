/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoleculeLookup.h" //Header spec.
#include "EnsemblePreprocessor.h" //For box total
#include "PRNG.h" //For index selection
#include "PDBSetup.h" //For init.
#include "Molecules.h" //For init.
#include <algorithm>
#include <utility>
#include <iostream>

void MoleculeLookup::Init(const Molecules& mols,
                          const pdb_setup::Atoms& atomData)
{
  numKinds = mols.GetKindsCount();
  molLookup = new uint[mols.count];
  molLookupCount = mols.count;

  //+1 to store end value
  boxAndKindStart = new uint[numKinds * BOX_TOTAL + 1];
  boxAndKindStartCount = numKinds * BOX_TOTAL + 1;

  // vector[box][kind] = list of mol indices for kind in box
  std::vector<std::vector<std::vector<uint> > > indexVector;
  indexVector.resize(BOX_TOTAL);
  fixedAtom.resize(mols.count);


  for (uint b = 0; b < BOX_TOTAL; ++b) {
    indexVector[b].resize(numKinds);
  }

  for(uint m = 0; m < mols.count; ++m) {
    uint box = atomData.box[atomData.startIdxRes[m]];
    uint kind = mols.kIndex[m];
    indexVector[box][kind].push_back(m);
    fixedAtom[m] = atomData.molBeta[m];

    //Find the kind that can be swap(beta == 0) or move(beta == 0 or 2)
    if(fixedAtom[m] == 0) {
      if(std::find(canSwapKind.begin(), canSwapKind.end(), kind) ==
          canSwapKind.end())
        canSwapKind.push_back(kind);

      if(std::find(canMoveKind.begin(), canMoveKind.end(), kind) ==
          canMoveKind.end())
        canMoveKind.push_back(kind);

    } else if(fixedAtom[m] == 2) {
      if(std::find(canMoveKind.begin(), canMoveKind.end(), kind) ==
          canMoveKind.end())
        canMoveKind.push_back(kind);

    }
  }


  uint* progress = molLookup;
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    for (uint k = 0; k < numKinds; ++k) {
      boxAndKindStart[b * numKinds + k] = progress - molLookup;
      progress = std::copy(indexVector[b][k].begin(),
                           indexVector[b][k].end(), progress);
    }
  }
  boxAndKindStart[numKinds * BOX_TOTAL] = mols.count;
}

uint MoleculeLookup::NumInBox(const uint box) const
{
  return boxAndKindStart[(box + 1) * numKinds]
         - boxAndKindStart[box * numKinds];
}

void MoleculeLookup::TotalAndDensity
(uint * numByBox, uint * numByKindBox, real * molFractionByKindBox,
 real * densityByKindBox, real const*const volInv) const
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    numByBox[b] = 0.0;
    for (uint k = 0; k < numKinds; ++k) {
      uint numMK = NumKindInBox(k, b);
      uint mkIdx = k + numKinds * b;
      numByKindBox[mkIdx] = numMK;
      densityByKindBox[mkIdx] = numMK * volInv[b];
      numByBox[b] += numMK;
    }

    //Calculate mol fractions
    if (numKinds > 1) {
      for (uint k = 0; k < numKinds; ++k) {
        uint mkIdx = k + numKinds * b;
        if (numByBox[b] > 0) {
          molFractionByKindBox[mkIdx] = (real)numByKindBox[mkIdx] /
                                        (real)numByBox[b];
        } else {
          molFractionByKindBox[mkIdx] = 0.0;
        }
      }
    }

  }

}

#ifdef VARIABLE_PARTICLE_NUMBER


bool MoleculeLookup::ShiftMolBox(const uint mol, const uint currentBox,
                                 const uint intoBox, const uint kind)
{
  uint index = std::find(
                 molLookup + boxAndKindStart[currentBox * numKinds + kind],
                 molLookup + boxAndKindStart[currentBox * numKinds + kind + 1], mol)
               - molLookup;
  assert(index != boxAndKindStart[currentBox * numKinds + kind + 1]);
  assert(molLookup[index] == mol);
  Shift(index, currentBox, intoBox, kind);
  return true;
}

void MoleculeLookup::Shift(const uint index, const uint currentBox,
                           const uint intoBox, const uint kind)
{
  uint oldIndex = index;
  uint newIndex;
  uint section = currentBox * numKinds + kind;
  if(currentBox >= intoBox) {
    while (section != intoBox * numKinds + kind) {
      newIndex = boxAndKindStart[section]++;
      uint temp = molLookup[oldIndex];
      molLookup[oldIndex] = molLookup[newIndex];
      molLookup[newIndex] = temp;
      oldIndex = newIndex;
      --section;
    }
  } else {
    while (section != intoBox * numKinds + kind) {
      newIndex = --boxAndKindStart[++section];
      uint temp = molLookup[oldIndex];
      molLookup[oldIndex] = molLookup[newIndex];
      molLookup[newIndex] = temp;
      oldIndex = newIndex;
    }
  }
}

#endif /*ifdef VARIABLE_PARTICLE_NUMBER*/


MoleculeLookup::box_iterator MoleculeLookup::box_iterator::operator++(int)
{
  box_iterator tmp = *this;
  ++(*this);
  return tmp;
}


MoleculeLookup::box_iterator::box_iterator(uint* _pLook, uint* _pSec)
  : pIt(_pLook + * _pSec) {}


MoleculeLookup::box_iterator MoleculeLookup::BoxBegin(const uint box) const
{
  return box_iterator(molLookup, boxAndKindStart + box * numKinds);
}

MoleculeLookup::box_iterator MoleculeLookup::BoxEnd(const uint box) const
{
  return box_iterator(molLookup, boxAndKindStart + (box + 1) * numKinds);
}
