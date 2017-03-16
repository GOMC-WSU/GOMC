#ifndef MOLECULELOOKUP_H
#define MOLECULELOOKUP_H

#include "EnsemblePreprocessor.h" //for BOX_TOTAL
#include "BasicTypes.h" //For uint

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

  MoleculeLookup() : molLookup(NULL), boxAndKindStart(NULL) {}

  ~MoleculeLookup()
  {
    delete[] molLookup;
    delete[] boxAndKindStart;
  }

  //Initialize this object to be consistent with Molecules mols
  void Init(Molecules const& mols, const pdb_setup::Atoms& atomData);

  uint GetNumKind(void) const
  {
    return numKinds;
  }

  //Returns number of given kind in given box
  uint NumKindInBox(const uint kind, const uint box) const;

  //!Returns total number of molecules in a given box
  uint NumInBox(const uint box) const;

  uint GetMolNum(const uint subIndex, const uint kind, const uint box)
  {
    return molLookup[boxAndKindStart[box * numKinds + kind]+subIndex];
  }

  void TotalAndDensity(uint * numByBox, uint * numByKindBox,
                       double * molFractionByBoxKind,
                       double * densityByKindBox,
		       double const * const volInv) const;

#ifdef VARIABLE_PARTICLE_NUMBER
  //Registers shift of mol into intoBox
  //Returns true if shift was successful, false otherwise
  bool ShiftMolBox(const uint mol, const uint currentBox,
                   const uint intoBox);

  //Same as above, but providing kind increases performance
  bool ShiftMolBox(const uint mol, const uint currentBox,
                   const uint intoBox, const uint kind);
#endif

  //iterator to traverse all the molecules in a particular box
  class box_iterator;
  friend class MoleculeLookup::box_iterator;
  box_iterator BoxBegin(const uint box) const;
  box_iterator BoxEnd(const uint box) const;

private:

#ifdef VARIABLE_PARTICLE_NUMBER
  void Shift(const uint index, const uint currentBox,
             const uint intoBox, const uint kind);
#endif


  //array of indices for type Molecule, sorted by box and kind for
  //move selection
  uint* molLookup;
  //index [BOX_TOTAL * kind + box] is the first element of that kind/box in
  //molLookup
  //index [BOX_TOTAL * kind + box + 1] is the element after the end
  //of that kind/box
  uint* boxAndKindStart;
  uint numKinds;
};

inline uint MoleculeLookup::NumKindInBox(const uint kind, const uint box) const
{
  return boxAndKindStart[box * numKinds + kind + 1] -
         boxAndKindStart[box * numKinds + kind];
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
  box_iterator(uint * _pLook, uint * _pSec);
  uint* pIt;
};

#endif
