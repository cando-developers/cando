/*
    File: energyNonbond.h
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
 *	energyNonbond.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyPeriodicBoundaryConditionNonbond_H  //[
#define	EnergyPeriodicBoundaryConditionNonbond_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energyNonbond.h>
//#include "geom/render.fwd.h"// energyNonbond.h wants DisplayList needs render.fwd.h

namespace chem {

num_real	_evaluateEnergyOnly_PeriodicBoundaryConditionNonbond(ScoringFunction_sp score,
                                                             num_real x1, num_real y1, num_real z1,
                                                             num_real x2, num_real y2, num_real z2,
                                                             num_real dA, num_real dC, num_real dQ1Q2 );


class EnergyPeriodicBoundaryConditionsNonbond_O : public EnergyNonbond_O
{
  LISP_CLASS(chem,ChemPkg,EnergyPeriodicBoundaryConditionsNonbond_O,"EnergyPeriodicBoundaryConditionsNonbond",EnergyNonbond_O);
  
  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m );
    
  virtual void evaluateUsingExcludedAtoms( ScoringFunction_sp score,
                                           NVector_sp 	pos,
                                           bool 		calcForce,
                                           gc::Nilable<NVector_sp> 	force,
                                           bool		calcDiagonalHessian,
                                           bool		calcOffDiagonalHessian,
                                           gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                           gc::Nilable<NVector_sp>	hdvec,
                                           gc::Nilable<NVector_sp> dvec);
  virtual void evaluateTerms( ScoringFunction_sp score,
                              NVector_sp 	pos,
                              bool 		calcForce,
                              gc::Nilable<NVector_sp> 	force,
                              bool		calcDiagonalHessian,
                              bool		calcOffDiagonalHessian,
                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                              gc::Nilable<NVector_sp>	hdvec,
                              gc::Nilable<NVector_sp> dvec);
  
  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(ScoringFunction_sp score, NVector_sp pos );

  // virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );
  virtual core::List_sp	checkForBeyondThresholdInteractionsWithPosition(ScoringFunction_sp score, NVector_sp pos, double threshold );

};

};

#endif //]
