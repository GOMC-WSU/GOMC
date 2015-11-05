/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef TRIALMOL_H
#define TRIALMOL_H

#include "../EnergyTypes.h"
#include "../XYZArray.h"
#include "../TransformMatrix.h"
#include "../../lib/BasicTypes.h"

class MoleculeKind;
class BoxDimensions;
class CalculateEnergy;

//! Class for keeping track of part-built molecules during CBMC
namespace cbmc
{

class TrialMol
{
   public:
      //!Construct TrialMol of kind k to be evaluated in box with axes ax.
      TrialMol(const MoleculeKind& k, const BoxDimensions& ax,
            uint box);
      //!Construct invalid default TrialMol
      TrialMol();

      TrialMol(const TrialMol& other);
      TrialMol& operator=(TrialMol other);
      friend void swap(TrialMol& a, TrialMol& b);

      //!True if this has been initialized to be valid
      bool IsValid() const { return (atomBuilt != NULL); }

      void AddAtom(uint index, const XYZ& position);

      void AddEnergy(const Energy& energy) { en += energy; }

      //!Confirms that atom at index i has been built (used for oldMols)
      void ConfirmOldAtom(uint i);

      //!Sets an orthonormal basis for coordinate conversion.
      /*!\param p1 Index of particle new additions will be bonded to
       * \param p2 Index of particle that will be in angles with new additions 
       * \param p3 Index of particle against which dihedrals will be measured
       */
      void SetBasis(uint p1, uint p2, uint p3);

      //!Sets an orthonormal basis for coordinate conversion.
      /*!\param p1 Index of particle new additions will be bonded to
       * \param p2 Index of particle that will be in angles with new additions 
       */
      void SetBasis(uint p1, uint p2);

      //!Shifts the current basis to the position of p1, but does not rotate it.
      void ShiftBasis(uint p1);

      //!Resets basis to box coordinate system
      void ResetBasis();

      //!Returns wrapped rectangular coordinates of a candidate position;
      XYZ GetRectCoords(double bond, double theta, double phi) const;

      XYZ RawRectCoords(double bond, double theta, double phi) const;
      
      // Returns the dihedral angle between two positions
      /* \param theta1 Theta spherical coordinate of first position
       * \param theta2 Theta spherical coordinate of second position
       * \param interior The interor angle between the positions
       */
      static double PhiBetweenAngles(double theta1, double theta2,
				     double interior);

      //!Return angle in radians between confirmed atoms a, b and c
      double GetTheta(uint a, uint b, uint c) const;

      //!Calculates theta and phi coords for atom in the current basis
      //!centered on lastAtom. theta in [0, pi], phi in (-pi, pi]
      void OldThetaAndPhi(uint atom, uint lastAtom, 
            double& theta, double& phi) const;

      const Energy& GetEnergy() const { return en; }
      double GetWeight() const { return totalWeight; }
      void SetWeight(double w) { totalWeight = w; }
      void MultWeight(double w) { totalWeight *= w; }

      uint GetBox() const { return box; }
      const BoxDimensions& GetAxes() const { return *axes; }
      const MoleculeKind& GetKind() const { return *kind; }

      //!Returns reference to coordinates of TrialMol.
      const XYZArray& GetCoords() const { return tCoords; }

      //!Returns position of atom i (undefined if it doesn't exist yet)
      XYZ AtomPosition(const uint atom) const { return tCoords.Get(atom); }

      //!Copies 1 molecule's worth of coordinates from coords[start] onwards
      void SetCoords(const XYZArray& coords, uint start);

      bool AtomExists(uint index) const { return atomBuilt[index]; }

      ~TrialMol();

   private:
      friend class CalculateEnergy;

      const MoleculeKind* kind;
      const BoxDimensions* axes;
      uint box;
      XYZArray tCoords;
      Energy en;
      double totalWeight;
      bool* atomBuilt;
      RotationMatrix growthToWorld;
      RotationMatrix worldToGrowth;
      XYZ basisPoint;
};
}




#endif

