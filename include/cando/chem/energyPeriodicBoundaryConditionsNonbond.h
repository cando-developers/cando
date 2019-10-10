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

#ifndef EnergyNonbond_H  //[
#define	EnergyNonbond_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
//#include "geom/render.fwd.h"// energyNonbond.h wants DisplayList needs render.fwd.h

namespace       chem {
namespace chem {

double	_evaluateEnergyOnly_PeriodicBoundaryConditionNonbond(double x1, double y1, double z1,
                                                             double x2, double y2, double z2,
                                                             double dA, double dC, double dQ1Q2 );


class EnergyPeriodicBoundaryConditionsNonbond_O : public EnergyNonbond_O
{
  LISP_CLASS(chem,ChemPkg,EnergyPeriodicBoundaryConditionsNonbond_O,"EnergyPeriodicBoundaryConditionsNonbond",EnergyNonbond_O);
  
  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m );
    
  virtual double evaluateAllComponent( ScoringFunction_sp scorer,
                                       NVector_sp 	pos,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec,
                                       gc::Nilable<NVector_sp> dvec);

  virtual void evaluateUsingExcludedAtoms( NVector_sp 	pos,
                                           bool 		calcForce,
                                           gc::Nilable<NVector_sp> 	force,
                                           bool		calcDiagonalHessian,
                                           bool		calcOffDiagonalHessian,
                                           gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                           gc::Nilable<NVector_sp>	hdvec,
                                           gc::Nilable<NVector_sp> dvec);
  virtual void evaluateTerms( NVector_sp 	pos,
                              bool 		calcForce,
                              gc::Nilable<NVector_sp> 	force,
                              bool		calcDiagonalHessian,
                              bool		calcOffDiagonalHessian,
                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                              gc::Nilable<NVector_sp>	hdvec,
                              gc::Nilable<NVector_sp> dvec);
  
  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos );

  // virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );
  virtual core::List_sp	checkForBeyondThresholdInteractionsWithPosition(NVector_sp pos, double threshold );

};

};

#endif //]
