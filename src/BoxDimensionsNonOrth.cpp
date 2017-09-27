/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.1
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "BoxDimensionsNonOrth.h"
#include "BoxDimensions.h"
#include "MoveConst.h" //For cutoff-related fail condition

void BoxDimensionsNonOrth::Init(config_setup::RestartSettings const& restart,
				config_setup::Volume const& confVolume,
				pdb_setup::Cryst1 const& cryst,
				double rc, double rcSq)
{
  rCut = rc;
  rCutSq = rcSq;
  for (uint b = 0; b < BOX_TOTAL; b++)
  {
    //if (restart.enable && cryst.hasVolume)
    //  axis = cryst.axis;
    if(confVolume.hasVolume)
      confVolume.axis[b].CopyRange(cellBasis[b], 0, 0, 3);
    else
    {
      fprintf(stderr,
	      "Error: Cell Basis not specified in PDB or in.dat files.\n");
      exit(EXIT_FAILURE);
    }
    //Find the length of a, b, c
    cellLength.Set(b, cellBasis[b].Length(0), cellBasis[b].Length(1),
		   cellBasis[b].Length(2));
    //Find Cosine Angle of alpha, beta and gamma
    cosAngle[b][0] = DotProduct(cellBasis[b].Get(1), cellBasis[b].Get(2)) /
      (cellLength.Get(b).y * cellLength.Get(b).z);
    cosAngle[b][1] = DotProduct(cellBasis[b].Get(0), cellBasis[b].Get(2)) /
      (cellLength.Get(b).x * cellLength.Get(b).z);
    cosAngle[b][2] = DotProduct(cellBasis[b].Get(0), cellBasis[b].Get(1)) /
      (cellLength.Get(b).x * cellLength.Get(b).y);
    //Calculate Cross Product
    XYZ axb = CrossProduct(cellBasis[b].Get(0), cellBasis[b].Get(1));
    XYZ bxc = CrossProduct(cellBasis[b].Get(1), cellBasis[b].Get(2));
    XYZ cxa = CrossProduct(cellBasis[b].Get(2), cellBasis[b].Get(0));
    //Calculate volume = A.(B x C)
    volume[b] = DotProduct(cellBasis[b].Get(0), bxc);
    volInv[b] = 1.0 / volume[b];
    //Calculate distance between two faces
    faceLength.Set(b, volume[b]/bxc.Length(), volume[b]/cxa.Length(),
		   volume[b]/axb.Length());
    //Calculate the adjoint and determinant
    double det = cellBasis[b].AdjointMatrix(cellBasis_Inv[b]);
    //Calculate the inverse matrix of cell basis
    cellBasis_Inv[b].ScaleRange(0, 3, 1.0/det);
    //Set the axis with unslant cell box
    //XYZ unslant = TransformUnSlant(cellLength.Get(b), b);
    //axis.Set(b, unslant.x, unslant.y, unslant.z); 
    axis.Set(b, cellLength[b]);
  }
  //We should consider the half of the face distance
  faceLength.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  halfAx.ScaleRange(0, BOX_TOTAL, 0.5);

  for (uint b = 0; b < BOX_TOTAL; b++)
  {
    //check to see if initial box size is cubic or not
    cubic[b] = ((axis.x[b] == axis.y[b]) && (axis.y[b] == axis.z[b]));
    //check to see if box is orthogonal or not
    orthogonal[b] = ((cosAngle[b][0] == 0.0) &&
		     (cosAngle[b][1] == 0.0) &&
		     (cosAngle[b][2] == 0.0));
  }
  constArea = confVolume.cstArea;
}

void BoxDimensionsNonOrth::CalcCellDimensions()
{
  for (uint b = 0; b < BOX_TOTAL; b++)
  {
    XYZ scale = axis.Get(b) / cellLength.Get(b);
    //Calculate new cell basis
    cellBasis[b].Scale(0, scale.x);
    cellBasis[b].Scale(1, scale.y);
    cellBasis[b].Scale(2, scale.z);
    //Set cell length
    cellLength.Set(b, axis[b]);
    //Calculate Cross Product
    XYZ axb = CrossProduct(cellBasis[b].Get(0), cellBasis[b].Get(1));
    XYZ bxc = CrossProduct(cellBasis[b].Get(1), cellBasis[b].Get(2));
    XYZ cxa = CrossProduct(cellBasis[b].Get(2), cellBasis[b].Get(0));
    //Calculate volume = A.(B x C)
    volume[b] = DotProduct(cellBasis[b].Get(0), bxc);
    volInv[b] = 1.0 / volume[b];
    //Calculate distance between two faces
    faceLength.Set(b, volume[b]/bxc.Length(), volume[b]/cxa.Length(),
		   volume[b]/axb.Length());
    //Calculate the adjoint and determinant
    double det = cellBasis[b].AdjointMatrix(cellBasis_Inv[b]);
    //Calculate the inverse matrix of cell basis
    cellBasis_Inv[b].ScaleRange(0, 3, 1.0/det);    
  }
  faceLength.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  halfAx.ScaleRange(0, BOX_TOTAL, 0.5);
}


BoxDimensionsNonOrth& BoxDimensionsNonOrth::operator=(BoxDimensionsNonOrth const& other)
{
  for (uint b = 0; b < BOX_TOTAL; ++b)
  {
    other.cellBasis[b].CopyRange(cellBasis[b], 0, 0, 3);
    other.cellBasis_Inv[b].CopyRange(cellBasis_Inv[b], 0, 0, 3);
    volume[b] = other.volume[b];
    volInv[b] = other.volInv[b];
    cubic[b] = other.cubic[b];
    orthogonal[b] = other.orthogonal[b];
    for(uint i = 0; i < 0; i++)
    {
      cosAngle[b][i] = other.cosAngle[b][i];
    }
  }
  other.axis.CopyRange(axis, 0, 0, BOX_TOTAL);
  other.halfAx.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  other.cellLength.CopyRange(cellLength, 0, 0, BOX_TOTAL);
  other.faceLength.CopyRange(faceLength, 0, 0, BOX_TOTAL);
  rCut = other.rCut;
  rCutSq = other.rCutSq;
  constArea = other.constArea;

  return *this;
}


void BoxDimensionsNonOrth::SetVolume(const uint b, const double vol)
{
   if (cubic[b] && !constArea)
   {
      double newAxX_b = pow(vol, (1.0/3.0));
      XYZ newAx_b(newAxX_b, newAxX_b, newAxX_b);
      axis.Set(b, newAx_b);
   }
   else if (constArea)
   {
     double area = axis.x[b] * axis.y[b];
     double newAxX_b = vol / area;
     XYZ newAx_b(axis.x[b], axis.y[b], newAxX_b);
     axis.Set(b, newAx_b);
   }
   else
   {
     double z_x = axis.z[b] / axis.x[b];
     double y_x = axis.y[b] / axis.x[b];
     double newAxX_b = pow(vol / (z_x * y_x), (1.0/3.0));
     XYZ newAx_b(newAxX_b, y_x * newAxX_b, newAxX_b * z_x);
     axis.Set(b, newAx_b);
   }
   //Calculate new cell dimension
   CalcCellDimensions();
}

XYZ BoxDimensionsNonOrth::MinImage(XYZ rawVec, const uint b) const
{
  rawVec.x = MinImageSigned(rawVec.x, faceLength.x[b], halfAx.x[b]);
  rawVec.y = MinImageSigned(rawVec.y, faceLength.y[b], halfAx.y[b]);
  rawVec.z = MinImageSigned(rawVec.z, faceLength.z[b], halfAx.z[b]);
  return rawVec;
}

void BoxDimensionsNonOrth::WrapPBC(double &x, double &y, double &z,
					  const uint b) const
{
  //convert XYZ to unslant
  XYZ unwrap(x, y, z);
  XYZ unslant = TransformUnSlant(unwrap, b);
  BoxDimensions::WrapPBC(unslant.x, axis.x[b]);
  BoxDimensions::WrapPBC(unslant.y, axis.y[b]);
  BoxDimensions::WrapPBC(unslant.z, axis.z[b]);
  //convert XYZ to slant
  XYZ slant = TransformSlant(unslant, b);
  x = slant.x;
  y = slant.y;
  z = slant.z; 
}

void BoxDimensionsNonOrth::UnwrapPBC(double & x, double & y, double & z,
					    const uint b, XYZ const& ref) const
{
  //convert XYZ to unslant
  XYZ wrap(x, y, z);
  XYZ unslant = TransformUnSlant(wrap, b);
  XYZ unslantRef = TransformUnSlant(ref, b);
  BoxDimensions::UnwrapPBC(unslant.x, unslantRef.x, axis.x[b], halfAx.x[b]);
  BoxDimensions::UnwrapPBC(unslant.y, unslantRef.y, axis.y[b], halfAx.y[b]);
  BoxDimensions::UnwrapPBC(unslant.z, unslantRef.z, axis.z[b], halfAx.z[b]);
  XYZ unwrap(x, y, z);
  //convert XYZ to slant
  XYZ slant = TransformSlant(unslant, b);
  x = slant.x;
  y = slant.y;
  z = slant.z; 
}

//Calculate transform
XYZ BoxDimensionsNonOrth::TransformUnSlant(const XYZ &A, const uint b) const
{
  XYZ temp;
  temp.x = A.x * cellBasis_Inv[b].Get(0).x + A.y * cellBasis_Inv[b].Get(1).x +
    A.z * cellBasis_Inv[b].Get(2).x;
  temp.y = A.x * cellBasis_Inv[b].Get(0).y + A.y * cellBasis_Inv[b].Get(1).y +
    A.z * cellBasis_Inv[b].Get(2).y;
  temp.z = A.x * cellBasis_Inv[b].Get(0).z + A.y * cellBasis_Inv[b].Get(1).z +
    A.z * cellBasis_Inv[b].Get(2).z;

  return temp;
}

//Calculate transform
XYZ BoxDimensionsNonOrth::TransformSlant(const XYZ &A,const uint b) const
{
  XYZ temp;
  temp.x = A.x * cellBasis[b].Get(0).x + A.y * cellBasis[b].Get(1).x +
    A.z * cellBasis[b].Get(2).x;
  temp.y = A.x * cellBasis[b].Get(0).y + A.y * cellBasis[b].Get(1).y +
    A.z * cellBasis[b].Get(2).y;
  temp.z = A.x * cellBasis[b].Get(0).z + A.y * cellBasis[b].Get(1).z +
    A.z * cellBasis[b].Get(2).z;

  return temp;
}

//Calculate AxB product
XYZ BoxDimensionsNonOrth::CrossProduct(const XYZ &A,const XYZ &B) const
{
  XYZ temp;
  temp.x = A.y * B.z - A.z * B.y;
  temp.y = A.z * B.x - A.x * B.z;
  temp.z = A.x * B.y - A.y * B.x;
  
  return temp;
}
