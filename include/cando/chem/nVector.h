/*
    File: nVector.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *      File:   nVector.h
 *
 *      Description:
 *              Routines to maintain N dimensional vectors
 */

#ifndef NVECTOR_H
#define NVECTOR_H

#include <clasp/core/common.h>
#include <cando/chem/nVector.fwd.h>
//#include "core/archiveNode.h"
#include <cando/chem/chemPackage.h>

namespace chem {
double	dotProduct( NVector_sp x, NVector_sp y, core::T_sp frozen );
  double	squared(NVector_sp me, core::T_sp frozen);
  double        magnitude(NVector_sp me, core::T_sp frozen);
  double	rmsMagnitude(NVector_sp me, core::T_sp frozen);
  double	angleWithVector(NVector_sp me, NVector_sp other, core::T_sp frozen);
  void          XPlusYTimesScalar( NVector_sp nvResult, NVector_sp nvX, NVector_sp nvD, double sigma, core::T_sp frozen);
  void	        copyVector(NVector_sp dest, NVector_sp orig);
  void	inPlaceAddTimesScalar( NVector_sp nvResult, NVector_sp nvY, double x, core::T_sp frozen);	// NVectorAdd(a,x,y)
  double	rmsDistanceFrom(NVector_sp u, NVector_sp v, core::T_sp frozen);
  NVector_sp  copy_nvector(NVector_sp original);
  
  void chem__nvector_ensure_identical(NVector_sp veca, NVector_sp vecb, double tol);
};

#endif          /* ifndef NVECTOR_H */
