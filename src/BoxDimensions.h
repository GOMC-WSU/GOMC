/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BOX_DIMENSIONS_H
#define BOX_DIMENSIONS_H

#include "EnsemblePreprocessor.h" //For BOX_TOTAL, ensembles
#include "../lib/BasicTypes.h" //For uint, etc.
#include "PDBSetup.h" //Primary source of volume.
#include "ConfigSetup.h" //Other potential source of volume (new sys only)
#include "XYZArray.h" //For axes
#include <cstdio>
#include <cstdlib>

#include <cassert>

//Use shortcuts when calculating Rcut
//#define RCUT_SHORTCUT


//TODO:
//Rename to PeriodicBoxes at some point.
struct BoxDimensions
{
 public:
   BoxDimensions() {}
   BoxDimensions(BoxDimensions const& other);
   BoxDimensions& operator=(BoxDimensions const& other);

   void Init(config_setup::RestartSettings const& restart,
	     config_setup::Volume const& confVolume, 
	     pdb_setup::Cryst1 const& cryst, double rc, double rcSq);

   XYZ GetAxis(const uint b) const { return axis.Get(b); }

   double GetTotVolume() const;

   void SetVolume(const uint b, const double vol);      
   
   uint ShiftVolume(BoxDimensions & newDim, double & scale, 
		    const uint b, const double delta) const;

   //!Calculate and execute volume exchange based on transfer
   /*!\param newDim return reference for new dimensions
    * \param scaleO return reference for linear scale factor of box bO
    * \param scaleN return reference for linear scale factor of box bN
    * \param bO box index for box O
    * \param bN box index for box N
    * \param transfer determines amount to be transfered
    */
   uint ExchangeVolume(BoxDimensions & newDim, double * scale, 
		       const double transfer) const;

   //Vector btwn two points, accounting for PBC, on an individual axis
   XYZ MinImage(XYZ rawVec, const uint b) const;
   
   //Wrap all coordinates in object.
   void WrapPBC(XYZArray & arr, const uint b) const;
   
   //Unwrap all coordinates in object.
   void UnwrapPBC(XYZArray & arr, const uint b, XYZ const& ref) const;
   
   //Wrap range of coordinates in object
   void WrapPBC(XYZArray & arr, const uint start, const uint stop,
		const uint b) const;

   //Unwrap range of coordinates in object
   void UnwrapPBC(XYZArray & arr, const uint start, const uint stop,
		  const uint b, XYZ const& ref) const;

   //Wrap one coordinate.
   XYZ WrapPBC(XYZ rawPos, const uint b) const;

   //Unwrap one coordinate.
   XYZ UnwrapPBC(XYZ& rawPos, const uint b, XYZ const& ref) const;
   
   //Unwrap one coordinate.
   void WrapPBC(double &x, double &y, double &z, const uint b) const
   { 
      WrapPBC(x, axis.x[b]);
      WrapPBC(y, axis.y[b]); 
      WrapPBC(z, axis.z[b]); 
   }

   //Unwrap one coordinate.
   void UnwrapPBC(double & x, double & y, double & z, 
		  const uint b, XYZ const& ref) const
   {
      UnwrapPBC(x, ref.x, axis.x[b], halfAx.x[b]); 
      UnwrapPBC(y, ref.y, axis.y[b], halfAx.y[b]); 
      UnwrapPBC(z, ref.z, axis.z[b], halfAx.z[b]); 
   }

   //Returns if within cutoff, if it is, gets distance -- 
   //with shortcut, same coordinate array
   bool InRcut(double & distSq, XYZ & dist, XYZArray const& arr, 
	       const uint i, const uint j, const uint b) const;
   
   //Dist squared -- with shortcut, two different coordinate arrays
   bool InRcut(double & distSq, XYZ & dist, XYZArray const& arr1, 
	       const uint i, XYZArray const& arr2, const uint j,
	       const uint b) const;

   //Returns if within cutoff, if it is, gets distance --
   //with shortcut, same coordinate array
   bool InRcut(double & distSq, XYZArray const& arr,
               const uint i, const uint j, const uint b) const;

   //Dist squared -- with shortcut, two different coordinate arrays
   bool InRcut(double & distSq, XYZArray const& arr1,
               const uint i, XYZArray const& arr2, const uint j,
               const uint b) const;

   bool InRcut(double distSq) const { return (distSq < rCutSq); }

   //Dist squared , two different coordinate arrays
   void GetDistSq(double & distSq, XYZArray const& arr1,
               const uint i, XYZArray const& arr2, const uint j,
               const uint b) const;
   //Dist squared with same coordinate array
   void GetDistSq(double & distSq, XYZArray const& arr, const uint i,
		  const uint j, const uint b) const;

   
   XYZArray axis;     //x, y, z dimensions of each box (a)
   XYZArray halfAx;   //x, y, z dimensions / 2 of each box (a)
   double volume[BOX_TOTAL];     //volume of each box in (a^3)
   double volInv[BOX_TOTAL];     //inverse volume of each box in (a^-3)

   double rCut;
   double rCutSq;
   double minBoxSize;
  
   //Dist. btwn two points, accounting for PBC, on an individual axis
   double MinImage(double& raw, const double ax, const double halfAx) const;
   double MinImageSigned(double raw, double ax, double halfAx) const;

   double WrapPBC(double& v, const double ax) const;

   double UnwrapPBC(double& v, const double ref,
		    const double ax, const double halfAx) const;
};

inline BoxDimensions::BoxDimensions(BoxDimensions const& other) : 
   axis(other.axis), halfAx(other.halfAx)
{
   for (uint b = 0; b < BOX_TOTAL; ++b)
   {
      volume[b] = other.volume[b];
      volInv[b] = other.volInv[b];
   }
   rCut = other.rCut;
   rCutSq = other.rCutSq;
}

inline BoxDimensions& BoxDimensions::operator=(BoxDimensions const& other)
{
   for (uint b = 0; b < BOX_TOTAL; ++b)
   {
      volume[b] = other.volume[b];
      volInv[b] = other.volInv[b];
   }
   other.axis.CopyRange(axis,0,0,BOX_TOTAL);
   other.halfAx.CopyRange(halfAx,0,0,BOX_TOTAL);
   minBoxSize = other.minBoxSize;
   rCut = other.rCut;
   rCutSq = other.rCutSq;
   return *this;
}

inline double BoxDimensions::GetTotVolume() const
{
   double sum = 0.0;
   for (uint b = 0; b < BOX_TOTAL; b++)
      sum += volume[b];
   return sum;
}


inline void BoxDimensions::SetVolume(const uint b, const double vol)
{
   double newAxX_b = pow(vol, (1.0/3.0));
   XYZ newAx_b(newAxX_b, newAxX_b, newAxX_b);
   volume[b] = vol;
   volInv[b] = 1.0/vol;
   axis.Set(b, newAx_b);
   newAx_b *= 0.5;
   halfAx.Set(b, newAx_b);
}

inline XYZ BoxDimensions::MinImage(XYZ rawVec, const uint b) const
{
   rawVec.x = MinImageSigned(rawVec.x, axis.x[b], halfAx.x[b]);
   rawVec.y = MinImageSigned(rawVec.y, axis.y[b], halfAx.y[b]);
   rawVec.z = MinImageSigned(rawVec.z, axis.z[b], halfAx.z[b]);
   return rawVec;
}

//Dist. btwn two points, accounting for PBC, on an individual axis
//
//Throws out sign (as per Brock's suggestion) as we don't care about it
//and thus can eliminate a branch and (potentially) one compare.
//
inline double BoxDimensions::MinImage
(double& raw, const double ax, const double halfAx) const
{
   raw = fabs(raw);
   //If shorter over periodic boundary, get that dist.
#ifdef NO_BRANCHING_MIN_IMAGE
   rawDiff = ax-raw;
   return (raw > halfAx)?rawDiff:raw;
#else
   if (raw > halfAx)
      raw = ax-raw;
   return raw; //...just pass back if distance is already minimum
#endif
}

inline double BoxDimensions::MinImageSigned(double raw, 
      double ax, double halfAx) const
{
   if (raw > halfAx)
      raw -= ax;
   else if (raw < -halfAx)
      raw += ax;
   return raw;
}

#ifdef RCUT_SHORTCUT
inline bool BoxDimensions::InRcut
(double & distSq, XYZ & dist,
 XYZArray const& arr, const uint i, const uint j, const uint b) const
{
   distSq = 0;
   dist.x = MinImageSigned(arr.DifferenceX(i, j), axis.x[b], halfAx.x[b]);
   if (rCut < dist.x)
      return false;
   dist.y = MinImageSigned(arr.DifferenceY(i, j), axis.y[b], halfAx.y[b]);
   if (rCut < dist.y)
      return false;
   distSq = dist.x*dist.x + dist.y*dist.y;
   if (rCutSq < distSq)
      return false;
   dist.z = MinImageSigned(arr.DifferenceZ(i, j), axis.z[b], halfAx.z[b]);
   if (rCut < dist.z)
      return false;
   distSq += dist.z*dist.z;
   return (rCutSq > distSq); 
}

inline bool BoxDimensions::InRcut
(double & distSq, XYZ & dist, XYZArray const& arr1, const uint i, 
 XYZArray const& arr2, const uint j, const uint b) const
{
   distSq = 0;
   dist.x = MinImageSigned(arr1.DifferenceX(i, arr2, j), 
			   axis.x[b], halfAx.x[b]);
   if (rCut < dist.x)
      return false;
   dist.y = MinImageSigned(arr1.DifferenceY(i, arr2, j), 
			   axis.y[b], halfAx.y[b]);
   if (rCut < dist.y)
      return false;
   distSq = dist.x*dist.x + dist.y*dist.y;
   if (rCutSq < distSq)
      return false;
   dist.z = MinImageSigned(arr1.DifferenceZ(i, arr2, j), 
			   axis.z[b], halfAx.z[b]);
   if (rCut < dist.z)
      return false;
   distSq += dist.z*dist.z;
   return (rCutSq > distSq); 
}

inline bool BoxDimensions::InRcut(double & distSq, XYZArray const& arr, 
				  const uint i, const uint j,
				  const uint b) const
{
   XYZ dist;
   distSq = 0;
   dist.x = MinImage(arr.DifferenceX(i, j), axis.x[b], halfAx.x[b]);
   if (rCut < dist.x)
      return false;
   dist.y = MinImage(arr.DifferenceY(i, j), axis.y[b], halfAx.y[b]);
   if (rCut < dist.y)
      return false;
   distSq = dist.x*dist.x + dist.y*dist.y;
   if (rCutSq < distSq)
      return false;
   dist.z = MinImage(arr.DifferenceZ(i, j), axis.z[b], halfAx.z[b]);
   if (rCut < dist.z)
      return false;
   distSq += dist.z*dist.z;
   return (rCutSq > distSq);
}

inline bool BoxDimensions::InRcut
(double & distSq, XYZArray const& arr1, const uint i,
 XYZArray const& arr2, const uint j, const uint b) const
{
   XYZ dist;
   distSq = 0;
   dist.x = MinImage(arr1.DifferenceX(i, arr2, j), axis.x[b], halfAx.x[b]);
   if (rCut < dist.x)
      return false;
   dist.y = MinImage(arr1.DifferenceY(i, arr2, j), axis.y[b], halfAx.y[b]);
   if (rCut < dist.y)
      return false;
   distSq = dist.x*dist.x + dist.y*dist.y;
   if (rCutSq < distSq)
      return false;
   dist.z = MinImage(arr1.DifferenceZ(i, arr2, j), axis.z[b], halfAx.z[b]);
   if (rCut < dist.z)
      return false;
   distSq += dist.z*dist.z;
   return (rCutSq > distSq);
}

#else

inline bool BoxDimensions::InRcut
(double & distSq, XYZ & dist,
 XYZArray const& arr, const uint i, const uint j, const uint b) const
{
   dist = MinImage(arr.Difference(i, j), b);
   distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
   return (rCutSq > distSq);
}


inline bool BoxDimensions::InRcut
(double & distSq, XYZ & dist, XYZArray const& arr1, const uint i,
 XYZArray const& arr2, const uint j, const uint b) const
{
   dist = MinImage(arr1.Difference(i, arr2, j), b);
   distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
   return (rCutSq > distSq);
}

inline bool BoxDimensions::InRcut
(double & distSq, XYZArray const& arr, const uint i, const uint j, const uint b) const
{
   XYZ dist = MinImage(arr.Difference(i, j), b);
   distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
   return (rCutSq > distSq);
}


inline bool BoxDimensions::InRcut
(double & distSq, XYZArray const& arr1, const uint i,
 XYZArray const& arr2, const uint j, const uint b) const
{
   XYZ dist = MinImage(arr1.Difference(i, arr2, j), b);
   distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
   return (rCutSq > distSq);
}


#endif

inline void BoxDimensions::GetDistSq(double & distSq, XYZArray const& arr1,
               const uint i, XYZArray const& arr2, const uint j,
               const uint b) const
{
   XYZ dist = MinImage(arr1.Difference(i, arr2, j), b);
   distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
}

inline void BoxDimensions::GetDistSq
(double & distSq, XYZArray const& arr, const uint i, const uint j,
 const uint b) const
{
   XYZ dist = MinImage(arr.Difference(i, j), b);
   distSq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
}

//Wrap one coordinate.
inline XYZ BoxDimensions::WrapPBC(XYZ rawPos, const uint b) const
{ 
   WrapPBC(rawPos.x, rawPos.y, rawPos.z, b);
   return rawPos;
}

//Wrap all coordinates in object.
inline void BoxDimensions::WrapPBC(XYZArray & arr, const uint b) const
{
   for (uint i = 0; i < arr.count; i++)
      WrapPBC(arr.x[i], arr.y[i], arr.z[i], b);
}

//Unwrap all coordinates in object.
inline void BoxDimensions::UnwrapPBC(XYZArray & arr, const uint b, XYZ 
      const& ref) const
{
   for (uint i = 0; i < arr.count; i++)
      UnwrapPBC(arr.x[i], arr.y[i], arr.z[i], b, ref);
}

//Wrap range of coordinates in object
inline void BoxDimensions::WrapPBC
(XYZArray & arr, const uint start, const uint stop, const uint b) const
{
   for (uint i = start; i < stop; i++)
      WrapPBC(arr.x[i], arr.y[i], arr.z[i], b);
}

//Unwrap range of coordinates in object
inline void BoxDimensions::UnwrapPBC(XYZArray & arr, const uint start, 
      const uint stop, const uint b,
      XYZ const& ref) const
{
   for (uint i = start; i < stop; i++)
      UnwrapPBC(arr.x[i], arr.y[i], arr.z[i], b, ref);
}

//
// Note, here we can't do the fabs trick as we need to know which end
// to wrap on.
//
inline double BoxDimensions::WrapPBC(double& v, const double ax) const
{
   //assert(v < 2*ax);

   //if ( v > ax ) //if +, wrap out to low end
   //   v -= ax;
   //else if ( v < 0 ) //if -, wrap to high end
   //   v += ax;
   //
   // Inspired by :
   // http://graphics.stanford.edu/~seander/bithacks.html#ConditionalNegate
   //
#ifdef NO_BRANCHING_WRAP
   // Inspired by :
   // http://graphics.stanford.edu/~seander/bithacks.html#ConditionalNegate
   //
   // Advantages:
   // No branching
   // 
   // Disadvantages:
   // Sometimes a couple of extra ops or a couple of extra compares.
   if (
   bool negate = (v > ax);
   double vNeg = v+(ax ^-negate)+negate;
   return (fabs(v-halfAx)>halfAx)?v:vNeg;
#else
   //Note: testing shows that it's most efficient to negate if true.
   //Source:
   // http://jacksondunstan.com/articles/2052
   if ( v >= ax ) //if +, wrap out to low end, on boundry will wrap to zero
      v -= ax;
   else if ( v < 0 ) //if -, wrap to high end
      v += ax;
   return v;
#endif
}

//Unwrap one coordinate.
inline XYZ BoxDimensions::UnwrapPBC(XYZ & rawPos, const uint b,
				    XYZ const& ref) const
{
   UnwrapPBC(rawPos.x, ref.x, axis.x[b], halfAx.x[b]); 
   UnwrapPBC(rawPos.y, ref.y, axis.y[b], halfAx.y[b]); 
   UnwrapPBC(rawPos.z, ref.z, axis.z[b], halfAx.z[b]); 
   return rawPos;
}

inline double BoxDimensions::UnwrapPBC
(double& v, const double ref, const double ax, const double halfAx) const
{
   //If absolute value of X dist btwn pt and ref is > 0.5 * box_axis
   //If ref > 0.5 * box_axis, add box_axis to pt (to wrap out + side)
   //If ref < 0.5 * box_axis, subtract box_axis (to wrap out - side)
   // uses bit hack to avoid branching for conditional
#ifdef NO_BRANCHING_UNWRAP
   bool negate = ( ref > halfAx );
   double vDiff = v+(ax^-negate)+negate;
   return (fabs(ref-v) > halfAx )?v:vDiff;
#else
   if (fabs(ref-v) > halfAx )
   {
      //Note: testing shows that it's most efficient to negate if true.
      //Source:
      // http://jacksondunstan.com/articles/2052
      if ( ref < halfAx )
	 v -= ax;
      else
	 v += ax;
   }
   return v;
#endif
}

#endif /*BOX_DIMENSIONS_H*/

