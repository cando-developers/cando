/*
    File: energyRigidBodyComponent.h
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
 *	energyRigidBodyComponent.h
 *
 */

#ifndef EnergyRigidBodyComponent_H  //[
#define	EnergyRigidBodyComponent_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>

namespace       chem {

  FORWARD(EnergyRigidBodyComponent);
class EnergyRigidBodyComponent_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyRigidBodyComponent_O,"EnergyRigidBodyComponent",EnergyComponent_O);
    DEFAULT_CTOR_DTOR(EnergyRigidBodyComponent_O);

  virtual	num_real evaluateAllComponent( ScoringFunction_sp scorer,
                                               NVector_sp 	pos,
                                               core::T_sp componentEnergy,
                                               bool 		calcForce,
                                               gc::Nilable<NVector_sp> 	force,
                                               bool		calcDiagonalHessian,
                                               bool		calcOffDiagonalHessian,
                                               gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                               gc::Nilable<NVector_sp>	hdvec,
                                               gc::Nilable<NVector_sp> dvec,
                                               core::T_sp activeAtomMask,
                                               core::T_sp debugInteractions ) = 0;

};







};

#endif //]
