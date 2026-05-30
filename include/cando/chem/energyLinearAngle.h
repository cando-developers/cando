/*
    File: energyLinearAngle.h
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
 *	energyLinearAngle.h
 *
 *	A one-sided linear-bend angle energy term for near-linear (sp) centers.
 *	E = kt*(1 + cos(theta)) (UFF linear bend) -- minimized at theta = 180 deg,
 *	harmonic in the deviation there, and singularity-free at linearity
 *	(written purely in cos(theta): no acos / no 1/sin(theta) factor), unlike
 *	the standard harmonic angle term in energyAngle.h.
 */

#ifndef EnergyLinearAngle_H  //[
#define	EnergyLinearAngle_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
#include <clasp/core/ql.h>
#include <cando/chem/energyKernels/linear_angle.h>


#include <cando/chem/chemPackage.h>


namespace chem {
class EnergyAtom;

inline	string	XmlTag_LinearAngle() { return "LinearAngle"; };
inline	string	XmlTag_EnergyLinearAngle() { return "EnergyLinearAngle"; };

/*! Store a linear-bend angle energy term.
 *
 */
class EnergyLinearAngle : public EnergyTerm
{
public:
	string	className()	{ return "EnergyLinearAngle"; };
public:
                // Parameters
	linear_angle_term	term;
                // Variables
        Atom_sp      	_Atom1;
        Atom_sp    	_Atom2;
        Atom_sp		_Atom3;
public:
  EnergyLinearAngle(Atom_sp a1, Atom_sp a2, Atom_sp a3, size_t i1, size_t i2, size_t i3, double kt) : term(kt,i1,i2,i3),_Atom1(a1),_Atom2(a2),_Atom3(a3) {};
  EnergyLinearAngle() {};

public:
        Atom_sp	getAtom1() { return this->_Atom1; };
        Atom_sp	getAtom2() { return this->_Atom2; };
        Atom_sp	getAtom3() { return this->_Atom3; };
	double	getKt()	{return this->term.kt;};
	double	getT();
	double	getTheta() {return this->getT();};
	double	getThetaDegrees(){return this->getT()/0.0174533;};

        core::List_sp encode() const;
        void decode(core::List_sp alist);
};
};

namespace translate {

template <>
struct	to_object<chem::EnergyLinearAngle >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::EnergyLinearAngle& angle)
  {
    return angle.encode();
  }
};

template <>
struct	from_object<chem::EnergyLinearAngle>
{
  typedef	chem::EnergyLinearAngle	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR("Implement me");
        }
};
};

namespace chem {

double	_evaluateEnergyOnly_LinearAngle(
    int I1, int I2, int I3, core::T_sp activeAtomMask,
    num_real x1, num_real y1, num_real z1,
    num_real x2, num_real y2, num_real z2,
    num_real x3, num_real y3, num_real z3,
    num_real kt );


FORWARD(EnergyLinearAngle);

class EnergyLinearAngle_O : public EnergyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyLinearAngle_O,"EnergyLinearAngle",EnergyComponent_O);
public:
  virtual bool restraintp() const override {return false;};
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public: // virtual functions inherited from Object
  void	initialize();
public:
  typedef EnergyLinearAngle	TermType;
public: // instance variables
  gctools::Vec0<TermType>	_Terms;
public:
  typedef gctools::Vec0<TermType>::iterator iterator;
  iterator begin() { return this->_Terms.begin(); };
  iterator end() { return this->_Terms.end(); };
  static EnergyLinearAngle_sp make(EnergyFunction_sp energyFunction);
public:
  virtual size_t numberOfTerms() { return this->_Terms.size();};
  CL_DEFMETHOD core::T_mv safe_amber_energy_linear_angle_term(size_t index) {
    if (index >= this->numberOfTerms() ) {
      SIMPLE_ERROR("Illegal term index {} must be less than {}" , index , this->_Terms.size() );
    }
    return Values(core::DoubleFloat_O::create(this->_Terms[index].term.kt),
                  core::make_fixnum(this->_Terms[index].term.i3x1),
                  core::make_fixnum(this->_Terms[index].term.i3x2),
                  core::make_fixnum(this->_Terms[index].term.i3x3)
                  );
  }

public:
  virtual std::string descriptionOfContents() const;
  void addTerm(const TermType& term);
  virtual void dumpTerms(core::HashTable_sp atomTypes);

  virtual core::List_sp extract_vectors_as_alist() const;
  virtual void fill_from_vectors_in_alist(core::List_sp vectors);


  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m,
                                          core::T_sp activeAtomMask );
  virtual double evaluateAllComponent( ScoringFunction_sp scorer,
                                       NVector_sp 	pos,
                                       core::T_sp energyScale,
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
      NVector_sp pos,
      core::T_sp activeAtomMask );

  void addLinearAngleTerm(AtomTable_sp at, Atom_sp a1, Atom_sp a2, Atom_sp a3, double kt);


  core::List_sp lookupLinearAngleTerms(AtomTable_sp at, Atom_sp a1, Atom_sp a2 , Atom_sp a3, core::HashTable_sp atomTypes );

  EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc);

  virtual void emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const;
  virtual size_t runTestCalls(core::T_sp stream, chem::NVector_sp pos) const;

  EnergyLinearAngle_O( const EnergyLinearAngle_O& ss ); //!< Copy constructor

  EnergyLinearAngle_O() : EnergyComponent_O() {};
};




};

#endif //]
