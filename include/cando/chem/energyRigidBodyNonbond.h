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


class RigidBodyNonbondCrossTerm {
public:
  REAL      dA;
  REAL      dC;
  RigidBodyNonbondCrossTerm(double a, double c) : dA(a), dC(c) {};
  RigidBodyNonbondCrossTerm() : dA(0.0), dC(0.0) {};
  core::List_sp encode() const;
  void decode(core::List_sp alist);
};

  
class RigidBodyAtomInfo {
public:
  core::T_sp   _Object;  // can be NIL
  INT          _TypeIndex;
  REAL         _Radius;
  REAL         _Epsilon;
  REAL         _Charge;
  Vector3      _Position;
  RigidBodyAtomInfo() : _Object(unbound<core::T_O>()) {};
  RigidBodyAtomInfo(core::T_sp a, double r, double e, double c, const Vector3& p) :
      _Object(a),
      _Radius(r),
      _Epsilon(e),
      _Charge(c),
      _Position(p) {};
  core::List_sp encode() const;
  void decode(core::List_sp alist);
};
  
SMART(EnergyRigidBodyNonbond);
};

namespace translate {

template <>
struct	to_object<chem::RigidBodyAtomInfo >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::RigidBodyAtomInfo& stretch)
  {
    return stretch.encode();
  }
};

template <>
struct	from_object<chem::RigidBodyAtomInfo>
{
  typedef	chem::RigidBodyAtomInfo	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR("Implement me");
        }
};


template <>
struct	to_object<chem::RigidBodyNonbondCrossTerm >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::RigidBodyNonbondCrossTerm& stretch)
  {
    return stretch.encode();
  }
};

template <>
struct	from_object<chem::RigidBodyNonbondCrossTerm>
{
  typedef	chem::RigidBodyNonbondCrossTerm	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR("Implement me");
        }
};
};


namespace chem {
class EnergyRigidBodyNonbond_O : public EnergyRigidBodyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyRigidBodyNonbond_O,"EnergyRigidBodyNonbond",EnergyRigidBodyComponent_O);
 public:
  CL_LISPIFY_NAME("make-energy-rigid-body-nonbond");
  CL_DEF_CLASS_METHOD static EnergyRigidBodyNonbond_sp make( core::Array_sp end_atoms );
public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
 public: // virtual functions inherited from Object
  void	initialize();
//	string	__repr__() const;
 public: // instance variables
  double		_DielectricConstant;
  double		_ScaleVdw;
  double		_ScaleElectrostatic;
  double		_EnergyVdw;
  double		_EnergyElectrostatic;
  core::SimpleVector_byte32_t_sp    _RigidBodyEndAtom;
  gctools::Vec0<RigidBodyAtomInfo>  _AtomInfoTable;
  size_t               _NumberOfTypes;
  gctools::Vec0<RigidBodyNonbondCrossTerm>  _CrossTerms;
 public:
  RigidBodyNonbondCrossTerm& crossTerm(size_t xi, size_t yi) { return this->_CrossTerms[this->_NumberOfTypes*xi+yi]; };
 public:
  void initializeCrossTerms(bool verbose);
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
  Vector3 getPosition(size_t index);
  CL_DEFMETHOD void energyRigidBodyNonbondSetTerm(gc::Fixnum index, core::T_sp object, double radius, double epsilon, double charge, const Vector3& position);
  
  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m,
                                          core::T_sp activeAtomMask );
    
  virtual double evaluateAllComponent( ScoringFunction_sp scorer,
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
                                         core::T_sp debugInteractions );


  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                               NVector_sp pos );

//    int countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn, core::LispPtr );

  virtual void dumpTerms(core::HashTable_sp atomTypes);
  core::List_sp parts_as_list(NVector_sp pos);
  size_t partsCoordinates(NVector_sp pos, size_t idx, core::SimpleVector_float_sp coords);
  
  core::ComplexVector_float_sp write_nonbond_atom_coordinates_to_complex_vector_float(NVector_sp pos, core::ComplexVector_float_sp parts);
  core::ComplexVector_float_sp write_rigid_body_coordinates_to_complex_vector_float(NVector_sp rigid_body_pos, core::Array_sp end_indicesx3, NVector_sp coordinates, core::ComplexVector_float_sp output);
  core::ComplexVector_sp write_nonbond_atoms_to_complex_vector(core::ComplexVector_sp parts);
  size_t number_of_nonbond_atoms() const;
 public:
 EnergyRigidBodyNonbond_O(core::SimpleVector_byte32_t_sp end_atoms) : _RigidBodyEndAtom(end_atoms) {
    this->resizeNonbondAtomInfoTable((*end_atoms)[end_atoms->length()-1]);
  }
};

};

#endif //]
