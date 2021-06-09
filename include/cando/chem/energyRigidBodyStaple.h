/*
    File: energyStretch.h
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
 *	energyRigidBodyStaple.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyRigidBodyStaple_H  //[
#define	EnergyRigidBodyStaple_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyRigidBodyComponent.h>

namespace       chem {


/*! Store a RigidBodyStaple energy term
 */
class EnergyRigidBodyStaple : public EnergyTerm {
public:
	string	className()	{ return "EnergyRigidBodyStaple"; };
public:
	REAL            ks;
        REAL            r0;
        INT             rigidBodyK;
        Vector3         pointK;
        INT             rigidBodyL;
        Vector3         pointL;
 EnergyRigidBodyStaple(double k, double r, size_t ck, const Vector3& pk, size_t cl, const Vector3& pl) :
        ks(k),
          r0(r),
          rigidBodyK(ck),
          pointK(pk),
          rigidBodyL(cl),
          pointL(pl) {};
  EnergyRigidBodyStaple() {};
         
};

 FORWARD(EnergyRigidBodyStaple);
 
class EnergyRigidBodyStaple_O : public EnergyRigidBodyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyRigidBodyStaple_O,"EnergyRigidBodyStaple",EnergyRigidBodyComponent_O);
 public:
 public: // virtual functions inherited from Object
  void	initialize();
  typedef EnergyRigidBodyStaple	TermType;
 public: // instance variables
  gctools::Vec0<TermType> _Terms;
 public:
  CL_LISPIFY_NAME("make-energy-rigid-body-staple");
  CL_DEF_CLASS_METHOD static EnergyRigidBodyStaple_sp make() {
    GC_ALLOCATE(EnergyRigidBodyStaple_O,energy);
    return energy;
  }
    
 public:
  CL_DEFMETHOD virtual size_t numberOfTerms() { return this->_Terms.size();};

  CL_DEFMETHOD void energy_rigid_body_staple_add_term(double ks, double r0, size_t cylk, const Vector3& pointk, size_t cyll, const Vector3& pointl);
    
  void addTerm(const TermType& term);
  virtual void dumpTerms();

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

  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos ) {IMPLEMENT_ME();};

//  virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos ) {IMPLEMENT_ME();};

  virtual string	beyondThresholdInteractionsAsString() {IMPLEMENT_ME();};

  core::List_sp parts_as_list(NVector_sp pos);

 public:
  EnergyRigidBodyStaple_O( const EnergyRigidBodyStaple_O& ss ); //!< Copy constructor

  DEFAULT_CTOR_DTOR(EnergyRigidBodyStaple_O);
};

};
#endif //]
