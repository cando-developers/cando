/*
    File: energyRigidBodyNonbond.h
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

#ifndef EnergyRigidBodyNonbond_H  //[
#define	EnergyRigidBodyNonbond_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/energyRigidBodyComponent.h>
//#include "geom/render.fwd.h"// energyNonbond.h wants DisplayList needs render.fwd.h

namespace       chem {


  struct RigidBodyNonbondCrossTerm {
    double      dA;
    double      dC;
  };
  
  struct RigidBodyAtomInfo {
    core::T_sp    _Atom;  // can be NIL
    size_t     _TypeIndex;
    double     _Radius;
    double     _Epsilon;
    double     _Charge;
    Vector3    _Position;
    RigidBodyAtomInfo() {};
  RigidBodyAtomInfo(core::T_sp a, double r, double e, double c, const Vector3& p) :
    _Atom(a),
      _Radius(r),
      _Epsilon(e),
      _Charge(c),
      _Position(p) {};
  };
  
SMART(EnergyRigidBodyNonbond);


#if 0
double	_evaluateEnergyOnly_Nonbond(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double dA, double dC, double dQ1Q2 );
#endif

class EnergyRigidBodyNonbond_O : public EnergyRigidBodyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyRigidBodyNonbond_O,"EnergyRigidBodyNonbond",EnergyRigidBodyComponent_O);
 public:
  CL_LISPIFY_NAME("make-energy-rigid-body-nonbond");
  CL_DEF_CLASS_METHOD static EnergyRigidBodyNonbond_sp make( core::Array_sp end_atoms );
 public: // virtual functions inherited from Object
  void	initialize();
//    void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;
 public: // instance variables
  double		_DielectricConstant;
  double		_ScaleVdw;
  double		_ScaleElectrostatic;
  double		_EnergyVdw;
  double		_EnergyElectrostatic;
  core::SimpleVector_byte32_t_sp    _RigidBodyEndAtom;
  gctools::Vec0<RigidBodyAtomInfo>  _AtomInfoTable;
 public:
  void	setDielectricConstant(double d) { this->_DielectricConstant = d; };
  double	getDielectricConstant() { return this->_DielectricConstant; };
  void	setVdwScale(double d) { this->_ScaleVdw = d; };
  double	getVdwScale()	{return this->_ScaleVdw; };
  void	setElectrostaticScale(double d) { this->_ScaleElectrostatic = d; };
  double	getElectrostaticScale()	{return this->_ScaleElectrostatic; };

  double	getVdwEnergy() { return this->_EnergyVdw; };
  double	getElectrostaticEnergy() { return this->_EnergyElectrostatic; };
  void resizeNonbondAtomInfoTable(size_t index) { this->_AtomInfoTable.resize(index); };

 public:
  void zeroEnergy();

  CL_DEFMETHOD void energyRigidBodyNonbondSetTerm(gc::Fixnum index, Atom_sp atom, double radius, double epsilon, double charge, const Vector3& position);
  
  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m );
    
  virtual double evaluateAll( NVector_sp 	pos,
                              bool 		calcForce,
                              gc::Nilable<NVector_sp> 	force,
                              bool		calcDiagonalHessian,
                              bool		calcOffDiagonalHessian,
                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                              gc::Nilable<NVector_sp>	hdvec,
                              gc::Nilable<NVector_sp> dvec);


  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                               NVector_sp pos );

//    int countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn, core::Lisp_sp );

  virtual	double	getEnergy();
  virtual void dumpTerms();
  core::List_sp parts_as_list(NVector_sp pos);
 public:
 EnergyRigidBodyNonbond_O(core::SimpleVector_byte32_t_sp end_atoms) : _RigidBodyEndAtom(end_atoms) {
    this->resizeNonbondAtomInfoTable((*end_atoms)[end_atoms->length()-1]);
  }
};

};

#endif //]
