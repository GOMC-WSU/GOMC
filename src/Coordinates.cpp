#include "Coordinates.h"
#include "TransformMatrix.h"
#include <algorithm>          //For copy
#include <cmath>
#include <cassert>

void Coordinates::InitFromPDB(pdb_setup::Atoms const& atoms)
{
  //Allocate master array and push stuff to it.
  XYZArray::Init(atoms.x.size());
  //Transfer without creating a new array.
  std::copy(atoms.x.begin(), atoms.x.end(), x);
  std::copy(atoms.y.begin(), atoms.y.end(), y);
  std::copy(atoms.z.begin(), atoms.z.end(), z);
  CheckCoordinate();
  comRef.CalcCOM();
}

void Coordinates::CheckCoordinate()
{
  int p, start, atom, length, stRange, endRange;
  XYZ min, max, diffV;

  for (uint b = 0; b < BOX_TOTAL; b++)
  {
    MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b),
      end = molLookRef.BoxEnd(b), endc = molLookRef.BoxEnd(b);
    //find the min and max coordinate
    stRange = molRef.MolStart(*thisMol);
    --endc;
    endRange = molRef.MolStart(*endc) + molRef.GetKind(*endc).NumAtoms();
    
    min.x = *std::min_element(x + stRange, x + endRange);
    max.x = *std::max_element(x + stRange, x + endRange);
    min.y = *std::min_element(y + stRange, y + endRange);
    max.y = *std::max_element(y + stRange, y + endRange);
    min.z = *std::min_element(z + stRange, z + endRange);
    max.z = *std::max_element(z + stRange, z + endRange);

    diffV = max - min;
    //printf("start: %d, end: %d\n", stRange, endRange);
    printf("Minimum coordinates in box %d: x= %4.3f, y= %4.3f, z= %4.3f\n",
	   b, min.x, min.y, min.z);
    printf("Maximum coordinates in box %d: x= %4.3f, y= %4.3f, z= %4.3f\n",
	   b, max.x, max.y, max.z);

    //check to see if molecules are in the box or not
    if( diffV.x > boxDimRef.axis.Get(b).x ||
	diffV.y > boxDimRef.axis.Get(b).y ||
	diffV.z > boxDimRef.axis.Get(b).z)
    {
      printf("Molecules are not packed inside the defined box dimension.\n");
      exit(0);
    }


    if(min.x < 0.0 || min.y < 0.0 || min.z < 0.0)
    {
      //shift all the molecules to positive axis
       XYZ shiftV = min;
       shiftV *= -1.0;
       printf("Note: Molecules in the box %d will be shifted to origin by \n vector [%4.3f, %4.3f, %4.3f].\n", b, shiftV.x, shiftV.y, shiftV.z);
    
       while (thisMol != end)
       {
	 start = molRef.MolStart(*thisMol);
	 MoleculeKind const& thisKind = molRef.GetKind(*thisMol);
      
	 for (p = 0; p < thisKind.NumAtoms(); p++)
	 {
	    atom = start + p;
	    x[atom] += shiftV.x;
	    y[atom] += shiftV.y;
	    z[atom] += shiftV.z;
	 }
	 ++thisMol;
       }
    
    }
  }

}

//Translate by a random amount
void Coordinates::TranslateRand
(XYZArray & dest, XYZ & newCOM,  uint & pStart, uint & pLen,
 const uint m, const uint b, const double max)
{
  XYZ shift = prngRef.SymXYZ(max);
  uint stop=0;
  //Get range.
  molRef.GetRange(pStart, stop, pLen, m);
  //Copy coordinates
  CopyRange(dest, pStart, 0, pLen);

  newCOM = comRef.Get(m);
  //Add translation
  dest.AddAll(shift);
  newCOM += shift;
  //Finish by rewrapping.
  boxDimRef.WrapPBC(dest, b);
  newCOM = boxDimRef.WrapPBC(newCOM, b);
}

//Rotate by a random amount.
void Coordinates::RotateRand
(XYZArray & dest, uint & pStart, uint & pLen,
 const uint m, const uint b, const double max)
{
  //Rotate (-max, max) radians about a uniformly random vector
  //Not uniformly random, but symmetrical wrt detailed balance
  RotationMatrix matrix = RotationMatrix::FromAxisAngle(
                            prngRef.Sym(max), prngRef.PickOnUnitSphere());

  XYZ center = comRef.Get(m);
  uint stop = 0;
  molRef.GetRange(pStart, stop, pLen, m);
  //Copy coordinates
  CopyRange(dest, pStart, 0, pLen);

  boxDimRef.UnwrapPBC(dest, b, center);
  //Do rotation
  for (uint p = 0; p < pLen; p++)   //Rotate each point.
  {
    dest.Add(p, -center);
    dest.Set(p, matrix.Apply(dest.Get(p)));
    dest.Add(p, center);
  }
  boxDimRef.WrapPBC(dest, b);

}

//scale all in each mol newCOM[m]/oldCOM[m]
void Coordinates::VolumeTransferTranslate
(uint & state, Coordinates & dest, COM & newCOM, BoxDimensions & newDim,
 COM const& oldCOM, const double max) const
{
  double scale[BOX_TOTAL], transfer = prngRef.Sym(max);
  for (uint b = 0; b < BOX_TOTAL; ++b)
  {
    scale[b] = 0.0;
  }
  //Scale cell
  state = boxDimRef.ExchangeVolume(newDim, scale, transfer);
  //If scaling succeeded (if it wouldn't take the box to below 2*rcut, cont.
  for (uint b = 0; b < BOX_TOTAL && (state == mv::fail_state::NO_FAIL); ++b)
  {
    TranslateOneBox(dest, newCOM, oldCOM, newDim, b, scale[b]);
  }
}



//Assumes dest is already initialized
void Coordinates::TranslateOneBox
(Coordinates & dest, COM & newCOM, COM const& oldCOM,
 BoxDimensions const& newDim, const uint b, const double scale) const
{
  uint pStart=0, pStop=0, pLen=0;
  MoleculeLookup::box_iterator curr = molLookRef.BoxBegin(b),
                               end = molLookRef.BoxEnd(b);
  while (curr != end)
  {
    molRef.GetRange(pStart, pStop, pLen, *curr);
    //Scale CoM for this molecule, translate all atoms by same amount
    newCOM.Scale(*curr, scale);
    XYZ shift = newCOM.Get(*curr);
    shift -= oldCOM.Get(*curr);
    //Translation of atoms in mol.
    //Unwrap coordinates
    XYZ oldCOMForUnwrap = oldCOM.Get(*curr);
    boxDimRef.UnwrapPBC(dest, pStart, pStop, b, oldCOMForUnwrap);
    dest.AddRange(pStart, pStop, shift);
    newDim.WrapPBC(dest, pStart, pStop, b);
    ++curr;
  }

}


