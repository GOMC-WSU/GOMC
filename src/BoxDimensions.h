/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BOX_DIMENSIONS_H
#define BOX_DIMENSIONS_H

#include "EnsemblePreprocessor.h" //For BOX_TOTAL, ensembles
#include "BasicTypes.h" //For uint, etc.
#include "PDBSetup.h" //Primary source of volume.
#include "ConfigSetup.h" //Other potential source of volume (new sys only)
#include "XYZArray.h" //For axes
#include "Forcefield.h"
#include "GeomLib.h"
#include <cstdio>
#include <cstdlib>

#include <cassert>

//Use shortcuts when calculating Rcut
//#define RCUT_SHORTCUT

class BoxDimensions
{
public:
  BoxDimensions()
  {
    for (uint b = 0; b < BOX_TOTAL; b++) {
      cellBasis[b] = XYZArray(3);
    }
  }
  BoxDimensions(BoxDimensions const& other);
  virtual BoxDimensions& operator=(BoxDimensions const& other);

  virtual void Init(config_setup::RestartSettings const& restart,
                    config_setup::Volume const& confVolume,
                    pdb_setup::Cryst1 const& cryst,
                    Forcefield const &ff);

  XYZ GetAxis(const uint b) const
  {
    XYZ temp;
    temp.x = axis[b][0];
    temp.y = axis[b][1];
    temp.z = axis[b][2];
    return temp;
  }

  double GetTotVolume(const uint b1, const uint b2) const;

  virtual void SetVolume(const uint b, const double vol);

  virtual uint ShiftVolume(BoxDimensions & newDim, XYZ & scale,
                           const uint b, const double delta) const;

  //!Calculate and execute volume exchange based on transfer
  virtual uint ExchangeVolume(BoxDimensions & newDim, XYZ * scale,
                              const double transfer, const uint *box) const;

  //Vector btwn two points, accounting for PBC, on an individual axis
  virtual XYZ MinImage(XYZ rawVec, const uint b) const;

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
  virtual void WrapPBC(double &x, double &y, double &z, const uint b) const;

  //Unwrap one coordinate.
  virtual void UnwrapPBC(double & x, double & y, double & z,
                         const uint b, XYZ const& ref) const;

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
/*
  bool InRcut(double distSq) const
  {
    return (distSq < rCutSq);
  }
  */

  //Dist squared , two different coordinate arrays
  void GetDistSq(double & distSq, XYZArray const& arr1,
                 const uint i, XYZArray const& arr2, const uint j,
                 const uint b) const;

  //Dist squared with same coordinate array
  void GetDistSq(double & distSq, XYZArray const& arr, const uint i,
                 const uint j, const uint b) const;
  
  //True if arr is inside cavDim with geometric center of center.
  bool InCavity(XYZ const& arr, XYZ const& center, XYZ const& cavDim,
                XYZArray const& invCav, const uint b) const;

  //Transform A to unslant coordinate
  virtual XYZ TransformUnSlant(const XYZ &A, const uint b) const;

  //Transform A to slant coordinate
  virtual XYZ TransformSlant(const XYZ &A, const uint b) const;

//private:
  //XYZArray axis;                  //x, y, z dimensions of each box (a)
  double axis[BOX_TOTAL][3];
  //XYZArray halfAx;                //x, y, z dimensions / 2 of each box (a)
  double halfAx[BOX_TOTAL][3];
  XYZArray cellBasis[BOX_TOTAL];  //x, y, z vector, 3 for each box
  double volume[BOX_TOTAL];       //volume of each box in (a^3)
  double volInv[BOX_TOTAL];       //inverse volume of each box in (a^-3)
  double cosAngle[BOX_TOTAL][3];  //alpha, beta, gamma for each box

  double rCut[BOX_TOTAL];
  double rCutSq[BOX_TOTAL];
  double minVol[BOX_TOTAL];

  bool cubic[BOX_TOTAL], orthogonal[BOX_TOTAL], constArea;

  //Dist. btwn two points, accounting for PBC, on an individual axis
  double MinImage(double& raw, const double ax, const double halfAx) const;
  double MinImageSigned(double raw, double ax, double halfAx) const;

  double WrapPBC(double& v, const double ax) const;

  double UnwrapPBC(double& v, const double ref,
                   const double ax, const double halfAx) const;
};


//Wrap one coordinate.
inline XYZ BoxDimensions::WrapPBC(XYZ rawPos, const uint b) const
{
  WrapPBC(rawPos.x, rawPos.y, rawPos.z, b);
  return rawPos;
}


//Unwrap one coordinate.
inline XYZ BoxDimensions::UnwrapPBC(XYZ & rawPos, const uint b,
                                    XYZ const& ref) const
{
  UnwrapPBC(rawPos.x, rawPos.y, rawPos.z, b, ref);
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
inline void BoxDimensions::WrapPBC(XYZArray & arr, const uint start,
                                   const uint stop, const uint b) const
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


inline bool BoxDimensions::InRcut(double & distSq, XYZ & dist,
                                  XYZArray const& arr, const uint i,
                                  const uint j, const uint b) const
{
  dist = MinImage(arr.Difference(i, j), b);
  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
  return (rCutSq[b] > distSq);
}


inline bool BoxDimensions::InRcut(double & distSq, XYZ & dist,
                                  XYZArray const& arr1, const uint i,
                                  XYZArray const& arr2, const uint j,
                                  const uint b) const
{
  dist = MinImage(arr1.Difference(i, arr2, j), b);
  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
  return (rCutSq[b] > distSq);
}

inline bool BoxDimensions::InRcut(double & distSq, XYZArray const& arr,
                                  const uint i, const uint j,
                                  const uint b) const
{
  XYZ dist = MinImage(arr.Difference(i, j), b);
  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
  return (rCutSq[b] > distSq);
}


inline bool BoxDimensions::InRcut(double & distSq, XYZArray const& arr1,
                                  const uint i, XYZArray const& arr2,
                                  const uint j, const uint b) const
{
  XYZ dist = MinImage(arr1.Difference(i, arr2, j), b);
  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
  return (rCutSq[b] > distSq);
}


inline void BoxDimensions::GetDistSq(double & distSq, XYZArray const& arr1,
                                     const uint i, XYZArray const& arr2,
                                     const uint j, const uint b) const
{
  XYZ dist = MinImage(arr1.Difference(i, arr2, j), b);
  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
}

inline void BoxDimensions::GetDistSq(double & distSq, XYZArray const& arr,
                                     const uint i, const uint j,
                                     const uint b) const
{
  XYZ dist = MinImage(arr.Difference(i, j), b);
  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
}

inline bool BoxDimensions::InCavity(XYZ const& arr, XYZ const& center,
                                    XYZ const& cavDim, XYZArray const& invCav,
                                    const uint b) const
{
  XYZ halfDim = cavDim * 0.5;
  halfDim *= halfDim;
  XYZ diff = MinImage(arr - center, b);
  diff = geom::Transform(invCav, diff);
  diff *= diff;
  if(diff.x > halfDim.x || diff.y > halfDim.y || diff.z > halfDim.z)
    return false;
  else
    return true;
}

//Calculate transform
inline XYZ BoxDimensions::TransformSlant(const XYZ &A, const uint b) const
{
  return A;
}

//Calculate transform
inline XYZ BoxDimensions::TransformUnSlant(const XYZ &A, const uint b) const
{
  return A;
}



#endif /*BOX_DIMENSIONS_H*/
