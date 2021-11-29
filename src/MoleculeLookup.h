/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOLECULELOOKUP_H
#define MOLECULELOOKUP_H

#include "EnsemblePreprocessor.h" //for BOX_TOTAL
#include "BasicTypes.h" //For uint
#include "algorithm"
#include "Forcefield.h"
#include <vector>
#include <cstring>
class CheckpointOutput;

namespace pdb_setup
{
class Atoms;
}
class Molecules;

//CLASS: MoleculeLookup
//This class provides a randomly chosen molecule index of a given kind in
//a given box requires updating with ShiftMolBox when a molecule is
//transferred from one box to another
class MoleculeLookup
{
public:

  MoleculeLookup(): molLookup(NULL), boxAndKindStart(NULL), boxAndKindSwappableCounts(NULL),
   molIndex(NULL), atomIndex(NULL), molKind(NULL), atomKind(NULL), atomCharge(NULL) {}


  ~MoleculeLookup()
  {
    delete[] molLookup;
    delete[] boxAndKindStart;
    delete[] molIndex;
    delete[] atomIndex;
    delete[] molKind;
    delete[] atomKind;
    delete[] atomCharge;
    delete[] boxAndKindSwappableCounts;
  }

  MoleculeLookup& operator=(const MoleculeLookup & rhs);
  bool operator==(const MoleculeLookup & rhs);


  static void swap(MoleculeLookup& oldMolLookup, MoleculeLookup& newMolLookup);


  //Initialize this object to be consistent with Molecules mols
  void Init(Molecules const& mols, const pdb_setup::Atoms& atomData,
            Forcefield &ff, bool restartFromCheckpoint);

  uint GetNumKind(void) const
  {
    return numKinds;
  }

  uint GetNumCanSwapKind(void) const
  {
    return canSwapKind.size();
  }

  uint GetNumCanMoveKind(void) const
  {
    return canMoveKind.size();
  }

  uint GetCanSwapKind(const uint k) const
  {
    return canSwapKind[k];
  }

  // returns true if molecule kind k can transfer between boxes
  bool IsKindSwapable(const uint k)
  {
    if(std::find(canSwapKind.begin(), canSwapKind.end(), k) == canSwapKind.end()) {
      return false;
    } else  {
      return true;
    }
  }

  uint GetCanMoveKind(const uint k) const
  {
    return canMoveKind[k];
  }

  // returns true if molecule kind k can move in it's box
  bool IsKindMoveable(const uint k)
  {
    if(std::find(canMoveKind.begin(), canMoveKind.end(), k) == canMoveKind.end()) {
      return false;
    } else  {
      return true;
    }
  }

  //Returns number of given kind in given box
  uint NumKindInBox(const uint kind, const uint box) const;

  uint NumKindInBoxSwappable(const uint kind, const uint box) const;

  //!Returns total number of molecules in a given box
  uint NumInBox(const uint box) const;

  uint GetBeta( const uint m) const
  {
    return fixedMolecule[m];
  }

  bool IsFix(const uint m) const
  {
    return (fixedMolecule[m] == 1);
  }

  bool IsNoSwap(const uint m) const
  {
    return (fixedMolecule[m] >= 1);
  }

  uint GetMolNum(const uint subIndex, const uint kind, const uint box)
  {
    return molLookup[boxAndKindStart[box * numKinds + kind] + subIndex];
  }

  uint GetSortedMolNum(const uint subIndex, const uint kind, const uint box)
  {
    return restartMoleculeIndices[molLookup[boxAndKindStart[box * numKinds + kind] + subIndex]];
  }

  // determine if molecule is in this box or not
  bool IsMoleculeInBox(const uint &molIdx, const uint &kindIdx, const uint &box) const
  {
    uint index = std::find(
                  molLookup + boxAndKindStart[box * numKinds + kindIdx],
                  molLookup + boxAndKindStart[box * numKinds + kindIdx + 1], molIdx)
                - molLookup;

    return ((molLookup[index] == molIdx));
  }

  // determine if atom is in this box or not// uses global atom index
  bool IsAtomInBox(const uint &aIdx, const uint &box)
  {
    return IsMoleculeInBox(molIndex[aIdx], molKind[aIdx], box);
  }

  void TotalAndDensity(uint * numByBox, uint * numByKindBox,
                       double * molFractionByBoxKind,
                       double * densityByKindBox,
                       double const * const volInv) const;

#ifdef VARIABLE_PARTICLE_NUMBER
  //Registers shift of mol into intoBox
  //Returns true if shift was successful, false otherwise
  bool ShiftMolBox(const uint mol, const uint currentBox,
                   const uint intoBox, const uint kind);
#endif

static uint GetConsensusMolBeta( const uint pStart,
                          const uint pEnd, 
                          const std::vector<double> & betas,
                          const uint m,
                          const uint box,
                          const std::string & name);

  //iterator to traverse all the molecules in a particular box
  class box_iterator;
  friend class MoleculeLookup::box_iterator;
  friend class CheckpointSetup;
  box_iterator BoxBegin(const uint box) const;
  box_iterator BoxEnd(const uint box) const;

//private:

#ifdef VARIABLE_PARTICLE_NUMBER
  void Shift(const uint index, const uint currentBox,
             const uint intoBox, const uint kind);
#endif

  // We keep the data stored in raw pointers for vectorization purposes.

  //array of indices for type Molecule, sorted by box and kind for
  //move selection
  uint32_t* molLookup;
  uint32_t molLookupCount;
  uint32_t atomCount;
  //index [BOX_TOTAL * kind + box] is the first element of that kind/box in
  //molLookup
  //index [BOX_TOTAL * kind + box + 1] is the element after the end
  //of that kind/box
  uint32_t* boxAndKindStart;
  uint32_t* boxAndKindSwappableCounts;
  uint32_t boxAndKindStartCount;
  uint32_t numKinds;
  int32_t *molIndex; // stores the molecule index for global atom index
  int32_t *atomIndex; // stores the local atom index for global atom index
  int32_t *molKind; // stores the molecule kind for global atom index
  int32_t *atomKind; // stores the atom kind for global atom index
  double *atomCharge; // stores the atom's charge for global

  // We save/load the data to/from vectors for serialization purposes.
  //array of indices for type Molecule, sorted by box and kind for
  //move selection
  std::vector<uint32_t> molLookupVec;
  //index [BOX_TOTAL * kind + box] is the first element of that kind/box in
  //molLookup
  //index [BOX_TOTAL * kind + box + 1] is the element after the end
  //of that kind/box
  std::vector<uint32_t> boxAndKindStartVec;
  std::vector<uint32_t> boxAndKindSwappableCountsVec;

  std::vector<int32_t> molIndexVec; // stores the molecule index for global atom index
  std::vector<int32_t> atomIndexVec; // stores the local atom index for global atom index
  std::vector<int32_t> molKindVec; // stores the molecule kind for global atom index
  std::vector<int32_t> atomKindVec; // stores the atom kind for global atom index
  std::vector<double> atomChargeVec; // stores the atom's charge for global atom index


  // For consistent trajectory ordering across checkpoints
  std::vector<uint32_t> restartMoleculeIndices;
  uint32_t restartedNumAtomsInBox[BOX_TOTAL];

  std::vector <uint32_t> fixedMolecule;
  std::vector <uint32_t> canSwapKind; //Kinds that can move intra and inter box
  std::vector <uint32_t> canMoveKind; //Kinds that can move intra box only

  // make CheckpointOutput class a friend so it can print all the private data
  friend class CheckpointOutput;
  private:
    friend class cereal::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & molLookupVec;
      ar & molLookupCount;
      ar & atomCount;
      ar & boxAndKindStartVec;
      ar & boxAndKindSwappableCountsVec;
      ar & boxAndKindStartCount;
      ar & numKinds;
      ar & restartMoleculeIndices;
      ar & restartedNumAtomsInBox;
      ar & fixedMolecule;
      ar & canSwapKind;
      ar & canMoveKind;
      ar & molIndexVec;
      ar & atomIndexVec;
      ar & molKindVec;
      ar & atomKindVec;
      ar & atomChargeVec;
    }


template <class T>
void wrapArrayInVector( T *sourceArray, size_t arraySize, std::vector<T, std::allocator<T> > &targetVector );
};

inline uint MoleculeLookup::NumKindInBox(const uint kind, const uint box) const
{
  return boxAndKindStart[box * numKinds + kind + 1] -
         boxAndKindStart[box * numKinds + kind];
}

inline uint MoleculeLookup::NumKindInBoxSwappable(const uint kind, const uint box) const
{
  return boxAndKindSwappableCounts[box * numKinds + kind];
}

class MoleculeLookup::box_iterator
{
  friend class MoleculeLookup;
public:
  bool operator== (const box_iterator& rhs) const
  {
    return (pIt == rhs.pIt);
  }
  bool operator!= (const box_iterator& rhs) const
  {
    return !(*this == rhs);
  }
  bool operator< (const box_iterator& rhs) const
  {
    return (pIt < rhs.pIt);
  }
  bool operator> (const box_iterator& rhs) const
  {
    return (rhs < *this);
  }
  bool operator<= (const box_iterator& rhs) const
  {
    return !(rhs < *this);
  }
  bool operator>= (const box_iterator& rhs) const
  {
    return !(*this < rhs);
  }
  box_iterator& operator++ ()
  {
    ++pIt;
    return *this;
  }
  box_iterator& operator-- ()
  {
    --pIt;
    return *this;
  }

  box_iterator operator++ (int);
  const uint& operator*() const
  {
    return *pIt;
  }
  box_iterator() : pIt(NULL) {}
private:
  box_iterator(const uint * _pLook, const uint * _pSec);
  uint const* pIt;
};


#endif
