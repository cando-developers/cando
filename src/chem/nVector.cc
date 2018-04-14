/*
    File: nVector.cc
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
       
#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/foundation.h>
#include <clasp/core/array.h>
#include <cando/chem/nVector.h>

namespace chem {
/*
 *      dotProduct
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the DOT product of two vectors.
 *      return:  (this) . x
 */
double	dotProduct( NVector_sp x, NVector_sp y )
{
  double          dDot;
  ASSERT(x->length()==y->length());
  dDot = 0.0;
  double* px = &(*x)[0];
  double* py = &(*y)[0];
  for ( size_t i(0), iEnd(x->length()); i<iEnd; ++i ) dDot += px[i]*py[i];
  return dDot;
}


double	squared(NVector_sp x)
{
  double          dDot;
  dDot = 0.0;
  double* dp = &(*x)[0];
  for ( size_t i(0), iEnd(x->length()); i<iEnd; ++i ) {
    dDot += dp[i]*dp[i];
  }
  return(dDot);
}

double magnitude(NVector_sp me)
{
  double  dDot;
  dDot = squared(me);
  return(sqrt(dDot));
}


double	rmsMagnitude(NVector_sp me)
{
  ASSERT(me->length()>0);
  double  dDot;
  dDot = squared(me);
  dDot /= (double)me->length();
  return(sqrt(dDot));
}


double	angleWithVector(NVector_sp me, NVector_sp other)
{
#define	VERY_SMALL 1.0e-6
  double	lenThis, lenOther, dot;
  lenThis = magnitude(me);
  lenOther = magnitude(other);
  dot = dotProduct(me->asSmartPtr(),other);
  if ( fabs(lenThis) > VERY_SMALL && fabs(lenOther)>VERY_SMALL) {
    dot /= (lenThis*lenOther);
    if ( dot > 1.0 ) dot = 1.0;
    if ( dot < -1.0 ) dot = -1.0;
    return acos(dot);
  }
  return 0.0;
}

void XPlusYTimesScalar( NVector_sp o, NVector_sp x, NVector_sp y, double s )
{
 LOG(BF("this->_Values.size() = %d") % o->length() );
 LOG(BF("x->size() = %d") % x->length() );
 LOG(BF("y->size() = %d") % y->length() );
 ASSERTP( o->length() == x->length(),"NVector_O::addTimesScalar>>mismatched length");
 ASSERT( o->length() == y->length());
 double* po = &(*o)[0];
 double* px = &(*x)[0];
 double* py = &(*y)[0];
 for ( int i(0),iEnd(o->length()); i<iEnd; ++i ) {
   po[i] = px[i]+py[i]*s;
 }
}

/*
 *      copyVector
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Copy the contents of nvVector into another.
 */
void	copyVector(NVector_sp dest, NVector_sp orig)
{
  ASSERTP(dest->length()==orig->length(),"NVector_O::copy>>mismatched length");
  double* ddest = &(*dest)[0];
  double* dorig = &(*orig)[0];
  for ( size_t i(0),iEnd(orig->length()); i<iEnd; ++i ) {
    ddest[i] = dorig[i];
  }
}

void inPlaceAddTimesScalar( NVector_sp result, NVector_sp dir, double s )
{_OF();
  ASSERT( result->length() == dir->length());
  double* presult = &(*result)[0];
  double* pdir = &(*dir)[0];
  for ( size_t i(0), iEnd(result->length()); i<iEnd; ++i ) {
    presult[i] += pdir[i]*s;
  }
}


double	rmsDistanceFrom(NVector_sp u, NVector_sp v)
{_OF();
  ASSERT(u->length() == v->length());
  double	e, sum;
  sum = 0.0;
  double* du = &(*u)[0];
  double* dv = &(*v)[0];
  for ( size_t i(0), iEnd(u->length()); i<iEnd; ++i ) {
    e = du[i]-dv[i];
    sum += e*e;
  }
  sum /= (float)(u->length());
  return sqrt(sum);
}

NVector_sp copy_nvector(NVector_sp orig )
{ 
  GC_COPY(NVector_O, vv, *orig );
  return vv;
}


};
