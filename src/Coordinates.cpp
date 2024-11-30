/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "Coordinates.h"

#include <cassert>

#include "TransformMatrix.h"

void Coordinates::InitFromPDB(pdb_setup::Atoms const &atoms) {
  // Allocate master array and push stuff to it.
  XYZArray::Init(atoms.beta.size());
  // Transfer without creating a new array.
  std::copy(atoms.x.begin(), atoms.x.end(), x);
  std::copy(atoms.y.begin(), atoms.y.end(), y);
  std::copy(atoms.z.begin(), atoms.z.end(), z);

  WrapCoordinate(atoms.min, atoms.max);
  comRef.CalcCOM();
}

void Coordinates::WrapCoordinate(const XYZ min[BOX_TOTAL],
                                 const XYZ max[BOX_TOTAL]) {
  int p, start, atom;

  bool sawZeroCoordinate;

  for (uint b = 0; b < BOX_TOTAL; b++) {
    sawZeroCoordinate = false;
    MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b),
                                 end = molLookRef.BoxEnd(b),
                                 endc = molLookRef.BoxEnd(b);
    /* Prevent segfault on empty boxes */
    if (thisMol == end) {
      printf("No molecules to wrap inside the simulation box %d:\n", b);
      continue;
    }

    printf("Minimum coordinates in box %d: x = %8.3f, y = %8.3f, z = %8.3f\n",
           b, min[b].x, min[b].y, min[b].z);
    printf("Maximum coordinates in box %d: x = %8.3f, y = %8.3f, z = %8.3f\n",
           b, max[b].x, max[b].y, max[b].z);

    printf("Wrapping molecules inside the simulation box %d:\n", b);
    while (thisMol != end) {
      start = molRef.MolStart(*thisMol);
      MoleculeKind const &thisKind = molRef.GetKind(*thisMol);

      for (p = 0; p < (int)thisKind.NumAtoms(); p++) {
        atom = start + p;
        if (!x[atom] && !y[atom] && !z[atom]) {
          if (sawZeroCoordinate) {
            printf("Error: Multiple atoms with zero coordinates were found.\n");
            exit(EXIT_FAILURE);
          } else {
            sawZeroCoordinate = true;
          }
        }

        boxDimRef.WrapPBC(x[atom], y[atom], z[atom], b);
        // check to see if it is in the box or not
        XYZ unSlant(x[atom], y[atom], z[atom]);
        unSlant = boxDimRef.TransformUnSlant(unSlant, b);

        if (unSlant.x > boxDimRef.axis.Get(b).x ||
            unSlant.y > boxDimRef.axis.Get(b).y ||
            unSlant.z > boxDimRef.axis.Get(b).z || unSlant.x < 0 ||
            unSlant.y < 0 || unSlant.z < 0) {
          printf(
              "Molecules %d is packed outside of the defined box dimension.\n",
              *thisMol);
          exit(EXIT_FAILURE);
        }
      }
      ++thisMol;
    }
  }
}

// Translate by a random amount
void Coordinates::TranslateRand(XYZArray &dest, XYZ &newCOM, uint &pStart,
                                uint &pLen, const uint m, const uint b,
                                const double max) {
  XYZ shift = prngRef.SymXYZ(max);
  uint stop = 0;
  // Get range.
  molRef.GetRange(pStart, stop, pLen, m);
  // Copy coordinates
  CopyRange(dest, pStart, 0, pLen);

  newCOM = comRef.Get(m);
  // Add translation
  dest.AddAll(shift);
  newCOM += shift;
  // Finish by rewrapping.
  boxDimRef.WrapPBC(dest, b);
  newCOM = boxDimRef.WrapPBC(newCOM, b);
}

// Rotate by a random amount
void Coordinates::RotateRand(XYZArray &dest, uint &pStart, uint &pLen,
                             const uint m, const uint b, const double max) {
  // Rotate (-max, max) radians about a uniformly random vector
  // Not uniformly random, but symmetrical wrt detailed balance
  // Note: In order for gcc to produce the same results as the Intel compiler,
  //       we need to create variables instead of using these calls to random
  //       functions as parameters to the FromAxisAngle function call.
  double theta = prngRef.Sym(max);
  XYZ axis = prngRef.PickOnUnitSphere();
  RotationMatrix matrix = RotationMatrix::FromAxisAngle(theta, axis);

  XYZ center = comRef.Get(m);
  uint stop = 0;
  molRef.GetRange(pStart, stop, pLen, m);
  // Copy coordinates
  CopyRange(dest, pStart, 0, pLen);

  boxDimRef.UnwrapPBC(dest, b, center);
  // Do rotation
  for (uint p = 0; p < pLen; p++) { // Rotate each point.
    dest.Add(p, -center);
    dest.Set(p, matrix.Apply(dest.Get(p)));
    dest.Add(p, center);
  }
  boxDimRef.WrapPBC(dest, b);
}

// scale all in each mol newCOM[m]/oldCOM[m]
void Coordinates::VolumeTransferTranslate(uint &state, Coordinates &dest,
                                          COM &newCOM, BoxDimensions &newDim,
                                          COM const &oldCOM, const double max,
                                          const uint *box) const {
  XYZ scale[2];
  double transfer = prngRef.Sym(max);

  // Scale cell
  state = boxDimRef.ExchangeVolume(newDim, scale, transfer, box);
  // If scaling succeeded (if it wouldn't take the box to below 2*rcut, cont.
  if (state == mv::fail_state::NO_FAIL) {
    for (uint b = 0; b < 2; ++b) {
      TranslateOneBox(dest, newCOM, oldCOM, newDim, b, scale[b]);
    }
  }
}

// Assumes dest is already initialized
void Coordinates::TranslateOneBox(Coordinates &dest, COM &newCOM,
                                  COM const &oldCOM,
                                  BoxDimensions const &newDim, const uint b,
                                  const XYZ &scale) const {
  uint pStart = 0, pStop = 0, pLen = 0;
  MoleculeLookup::box_iterator curr = molLookRef.BoxBegin(b),
                               end = molLookRef.BoxEnd(b);
  XYZ shift, oldCOMForUnwrap;
  XYZ unslant, slant;

  while (curr != end) {
    molRef.GetRange(pStart, pStop, pLen, *curr);
    // Scale CoM for this molecule, translate all atoms by same amount
    // convert the COM to unslant coordinate
    unslant = boxDimRef.TransformUnSlant(newCOM.Get(*curr), b);
    // scale the COM
    unslant *= scale;
    // convert to slant coordinate
    slant = newDim.TransformSlant(unslant, b);
    // calculate the difference of new and old COM
    newCOM.Set(*curr, slant);
    shift = newCOM.Get(*curr);
    shift -= oldCOM.Get(*curr);
    // Translation of atoms in mol.
    // Unwrap coordinates
    oldCOMForUnwrap = oldCOM.Get(*curr);
    boxDimRef.UnwrapPBC(dest, pStart, pStop, b, oldCOMForUnwrap);
    dest.AddRange(pStart, pStop, shift);
    newDim.WrapPBC(dest, pStart, pStop, b);
    ++curr;
  }
}
