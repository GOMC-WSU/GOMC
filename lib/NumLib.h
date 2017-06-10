#ifndef NUMERIC_LIB_H
#define NUMERIC_LIB_H

#include <limits> //for double limits
#include <vector> //for vector average
#include "BasicTypes.h" //For uint, XYZ

#define DBL_MAX 1.7976931348623158e+308

namespace num
{
   static const double dbl_margin = 0.00001;
   static const double qqFact = 167000.00;
   static const double BIGNUM = DBL_MAX;
   static const uint VDW_STD_KIND = 0, VDW_SHIFT_KIND = 1, VDW_SWITCH_KIND = 2;

   template <typename T>
   inline void BoundGt(double & val, const double bound)
   { if (val > bound) val = bound; }

   template <typename T>
   inline void BoundLt(double & val, const double bound)
   { if (val < bound) val = bound; }

   template <typename T>
   inline void BoundNZDecimal(T & val, const int mult)
   { BoundLt<T>(val,std::numeric_limits<T>::min() * mult); }

   template <typename T>
   inline void Bound(T & val, const double lower, const double upper)
   {
      BoundLt<T>(val, lower);
      BoundGt<T>(val, upper);
   }

   //Arithmetic mean.
   inline double MeanA(const uint v1, const uint v2) { return (v1+v2)/2.0; }
   //Geometric mean.
   inline double MeanG(const double v1, const double v2) { return sqrt(v1*v2); }
   //Arithmetic mean.
   inline double MeanA(std::vector<double> const& v1,
		       std::vector<double> const& v2,
                       const uint ix1, const uint ix2)
   { return (v1[ix1]+v2[ix2])*0.5; }
   //Arithmetic mean.
   inline double MeanA(std::vector<uint> const& v1,
                       std::vector<uint> const& v2,
                       const uint ix1, const uint ix2)
   {
#ifdef MIE_INT_ONLY
     return (v1[ix1]+v2[ix2])/2;
#else
     return ((double)(v1[ix1]+v2[ix2]))/2.0;
#endif
   }
   //Geometric mean.
   inline double MeanG(std::vector<double> const& v1,
		       std::vector<double> const& v2,
                       const uint ix1, const uint ix2)
   { return sqrt(v1[ix1]*v2[ix2]); }

   template <class Type>
   inline Type Sq(const Type v) { return v * v; }
   template <class Type>
   inline Type Cb(const Type v) { return v * v * v; }
   template <class Type>
   inline void Cb(Type & s, Type & c, const Type v)
   { s = v * v; c = s * v; }

   inline double POW(const double d2, const double d4, const double d6,
		     uint e)
   {
      double result = (e & 0x1 ? sqrt(d2) : 1.0);
      e >>= 1;
      switch (e)
      {
	 case 0: break;
	 case 1: result *= d2; break;
	 case 2: result *= d4; break;
	 case 3: result *= d6; break;
	 case 4: result *= Sq(d4); break;
	 case 5: result *= d4 * d6; break;
	 case 6: result *= Sq(d6); break;
	 case 7: result *= Sq(d6) * d2; break;
	 case 8: result *= Sq(d6) * d4; break;
	 case 9: result *= Sq(d6) * d6; break;
	 case 10: result *= Sq(d6) * Sq(d4); break;
	 case 11: result *= Sq(Sq(d4)) * d6; break;
	 case 12: result *= Sq(Sq(d6)); break;
	 case 13: result *= Sq(Sq(d6))*d2; break;
	 case 14: result *= Sq(Sq(d6))*d4; break;
	 case 15: result *= Sq(Sq(d6))*d6; break;
	 case 16: result *= Sq(Sq(Sq(d4))); break;
	 case 17: result *= Sq(Sq(Sq(d4)))*d2; break;
	 case 18: result *= Sq(Sq(d6)*d6); break;
	 case 19: result *= Sq(Sq(d6)*d6)*d2; break;
	 case 20: result *= Sq(Sq(d6*d4)); break;
	 case 21: result *= Sq(Sq(d6*d4))*d2; break;
	 case 22: result *= Sq(Sq(Sq(d4)) * d6); break;
	 case 23: result *= Sq(Sq(Sq(d4)) * d6)*d2; break;
	 case 24: result *= Sq(Sq(Sq(d6))); break;
	 case 25: result *= Sq(Sq(Sq(d6)))*d2; break;
      }
      return result;
   }
}

#endif /*NUMERIC_LIB_H*/
