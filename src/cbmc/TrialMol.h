/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef TRIALMOL_H
#define TRIALMOL_H

#include "BasicTypes.h"
#include "EnergyTypes.h"
#include "GeomLib.h"
#include "MoleculeKind.h"
#include "TransformMatrix.h"
#include "XYZArray.h"

class MoleculeKind;
class BoxDimensions;
class CalculateEnergy;

//! Class for keeping track of part-built molecules during CBMC
namespace cbmc {
struct Bonds {
public:
  Bonds(const BondList &bList) {
    for (uint b = 0; b < bList.count; b++) {
      a0.push_back(bList.part1[b]);
      a1.push_back(bList.part2[b]);
      kind.push_back(bList.kinds[b]);
      bondExist.push_back(false);
      count = bList.count;
    }
  }

  Bonds() { count = 0; }

  Bonds(const Bonds &other) {
    a0 = other.a0;
    a1 = other.a1;
    kind = other.kind;
    bondExist = other.bondExist;
    count = other.count;
  }

  void AddBond(const uint p0, const uint p1) {
    for (uint b = 0; b < count; b++) {
      if (a0[b] == p0 && a1[b] == p1) {
        bondExist[b] = true;
      } else if (a0[b] == p1 && a1[b] == p0) {
        bondExist[b] = true;
      }
    }
  }

  bool BondExist(const uint p0, const uint p1) {
    for (uint b = 0; b < count; b++) {
      if (a0[b] == p0 && a1[b] == p1) {
        return bondExist[b];
      } else if (a0[b] == p1 && a1[b] == p0) {
        return bondExist[b];
      }
    }
    // We should not reach here
    return false;
  }

  void Unset() { std::fill(bondExist.begin(), bondExist.end(), false); }

private:
  std::vector<uint> a0, a1, kind;
  std::vector<bool> bondExist;
  uint count;
};

class TrialMol {
public:
  //! Construct TrialMol of kind k to be evaluated in box with axes ax.
  TrialMol(const MoleculeKind &k, const BoxDimensions &ax, uint box);
  //! Construct invalid default TrialMol
  TrialMol();

  TrialMol(const TrialMol &other);
  TrialMol &operator=(TrialMol other);
  friend void swap(TrialMol &a, TrialMol &b);

  //! True if this has been initialized to be valid
  bool IsValid() const { return (atomBuilt != NULL); }

  void AddAtom(uint index, const XYZ &position);
  void SetAtomCoords(uint index, const XYZ &loc);

  void AddEnergy(const Energy &energy) { en += energy; }

  // Keep the tcoords and reset everythings
  void Reset() {
    totalWeight = 1.0;
    std::fill_n(atomBuilt, kind->NumAtoms(), false);
    growthToWorld.LoadIdentity();
    en.Zero();
    bonds.Unset();
    overlap = false;
  }

  //! Confirms that atom at index i has been built (used for oldMols)
  void ConfirmOldAtom(uint i);

  //! Sets an orthonormal basis for coordinate conversion.
  /*!\param p1 Index of particle new additions will be bonded to
   * \param p2 Index of particle that will be in angles with new additions
   * \param p3 Index of particle against which dihedrals will be measured
   */
  void SetBasis(uint p1, uint p2, uint p3);

  //! Sets an orthonormal basis for coordinate conversion.
  /*!\param p1 Index of particle new additions will be bonded to
   * \param p2 Index of particle that will be in angles with new additions
   */
  void SetBasis(uint p1, uint p2);

  //! Shifts the current basis to the position of p1, but does not rotate it.
  void ShiftBasis(uint p1);

  //! Shifts the current basis to the XYZ coordinate.
  void ShiftBasis(XYZ cent);

  //! Resets basis to box coordinate system
  void ResetBasis();

  //! Returns wrapped rectangular coordinates of a candidate position;
  XYZ GetRectCoords(double bond, double theta, double phi) const;

  XYZ RawRectCoords(double bond, double theta, double phi) const;

  // Returns the dihedral angle between two positions
  /* \param theta1 Theta spherical coordinate of first position
   * \param theta2 Theta spherical coordinate of second position
   * \param interior The interor angle between the positions
   */
  static double PhiBetweenAngles(double theta1, double theta2, double interior);

  //! Return angle in radians between confirmed atoms a-b-c
  double GetTheta(uint a, uint b, uint c) const;

  //! Return dihedral in radians between confirmed atoms a-b-c-d
  double GetPhi(uint a, uint b, uint c, uint d) const;

  //! Calculates theta and phi coords for atom in the current basis
  //! centered on lastAtom. theta in [0, pi], phi in (-pi, pi]
  void OldThetaAndPhi(uint atom, uint lastAtom, double &theta,
                      double &phi) const;

  //! calculate distance between atoms belong to specified angle
  double AngleDist(const double b1, const double b2, const double theta);

  //! calculate distance between atoms belong to specified dihedral
  double DihedDist(const double b1, const double b2, const double b3,
                   const double theta1, const double theta2, const double phi);

  //! calculate distance between two atom in oldMol
  double OldDistSq(const uint atom, const uint lastAtom);

  // calculate min image distance between a and b
  double DistSq(const XYZ &a, const XYZ &b);

  const Energy &GetEnergy() const { return en; }
  double GetWeight() const { return totalWeight; }
  void SetWeight(double w) { totalWeight = w; }
  void MultWeight(double w) { totalWeight *= w; }

  uint GetBox() const { return box; }
  const BoxDimensions &GetAxes() const { return *axes; }
  const MoleculeKind &GetKind() const { return *kind; }

  bool OneFour() const { return kind->oneFour; }

  //! Returns reference to coordinates of TrialMol.
  const XYZArray &GetCoords() const { return tCoords; }

  //! Returns reference to coordinates of the backup TrialMol.
  const XYZArray &GetBCoords() const { return bCoords; }

  //! Returns position of atom i (undefined if it doesn't exist yet)
  XYZ AtomPosition(const uint atom) const { return tCoords.Get(atom); }

  //! Copies 1 molecule's worth of coordinates from coords[start] onwards to
  //! tCoords
  void SetCoords(const XYZArray &coords, uint start);

  //! Copies 1 molecule's worth of coordinates from coords[start] onwards to
  //! bCoords
  void SetBCoords(const XYZArray &coords, uint start);

  bool AtomExists(uint index) const { return atomBuilt[index]; }

  void UpdateOverlap(const bool state) { overlap |= state; }

  bool HasOverlap() const { return overlap; }

  // Used in MEMC move
  void SetSeed(const XYZ &coords, const XYZ &cav, const bool inCav,
               const bool fixCOM, const bool rotBB);
  void SetSeed(const bool inCav, const bool fixCOM, const bool rotBB);
  void SetBackBone(const int bb[2]);
  // sets atom index where molecule start growing
  void SetGrowingAtomIndex(const int &idx) { growingAtomIndex = idx; }
  XYZ Transform(const XYZ &a) { return geom::Transform(cavMatrix, a); }
  void TransposeMatrix(XYZArray &invMatrix) {
    return geom::TransposeMatrix(invMatrix, cavMatrix);
  }
  bool HasCav() const { return comInCav; }
  bool COMFix() const { return comFix; }
  bool RotateBB() const { return rotateBB; }
  void SetCavMatrix(const XYZArray &matrix);
  XYZ GetCavityCenter() const { return cavityCenter; }
  XYZ GetCavity() const { return cavity; }
  // return unwrap com of tcoords so tcoords must be set
  XYZ GetCOM();

  // returns the backbone index
  uint GetAtomBB(const uint i) const { return backbone[i]; }

  // returns the growing atom index
  int GetGrowingAtomIndex() const { return growingAtomIndex; }

  // set built bond to true
  void AddBonds(const uint p0, const uint p1) { bonds.AddBond(p0, p1); }
  // check to see if bond exist or not
  bool BondsExist(const uint p0, const uint p1) {
    return bonds.BondExist(p0, p1);
  }

  ~TrialMol();

private:
  friend class CalculateEnergy;

  const MoleculeKind *kind;
  const BoxDimensions *axes;
  uint box;
  XYZArray tCoords, cavMatrix;
  XYZArray bCoords; // used to find the angle and theta in rings molecule
  Energy en;
  double totalWeight;
  RotationMatrix growthToWorld;
  RotationMatrix worldToGrowth;
  XYZ basisPoint;
  XYZ cavityCenter, cavity; // The center and cavity dimensions
  int backbone[2];
  int growingAtomIndex; // use to start growing atom using CD-CBMC
  bool comInCav, comFix, rotateBB;
  bool overlap;
  bool *atomBuilt;
  // To check the status of built bonds
  Bonds bonds;
};
} // namespace cbmc

#endif
