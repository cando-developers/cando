/*
    File: energyAngle.h
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
 *	energyAngle.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyAngle_H  //[
#define	EnergyAngle_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>



#include <cando/chem/chemPackage.h>


namespace chem {
class EnergyAtom;

SMART(FFAngle);

inline	string	XmlTag_Angle() { return "Angle"; };
inline	string	XmlTag_EnergyAngle() { return "EnergyAngle"; };

struct TermAngle {
	REAL		kt;	//!< Angle bend force constant, this must match Mathematica code!
	REAL		t0;	//!< Angle bend equilibrium angle (Radians), this must match Mathematica code!
	INT		I1;	//!< i*3 index into coordinate vector for atom1, this must match Mathematica code!
	INT		I2;	//!< i*3 index into coordinate vector for atom2, this must match Mathematica code!
	INT		I3;	//!< i*3 index into coordinate vector for atom3, this must match Mathematica code!
  TermAngle(REAL k, REAL t, INT i1, INT i2, INT i3) : kt(k), t0(t), I1(i1), I2(i2), I3(i3) {};
  TermAngle() {};
};

/*! Store an Angle energy term.
 *
 */
class EnergyAngle : public EnergyTerm 
{
public:
	string	className()	{ return "EnergyAngle"; };
public:
                // Parameters
	TermAngle	term;
                // Variables
        Atom_sp      	_Atom1;
        Atom_sp    	_Atom2;
        Atom_sp		_Atom3;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_Angle_debugEvalDeclares.cc>
#endif
public:
  EnergyAngle(Atom_sp a1, Atom_sp a2, Atom_sp a3, size_t i1, size_t i2, size_t i3, double kt, double t0) : term(kt,t0,i1,i2,i3),_Atom1(a1),_Atom2(a2),_Atom3(a3) {};
  EnergyAngle() {};

public:
        Atom_sp	getAtom1() { return this->_Atom1; };
        Atom_sp	getAtom2() { return this->_Atom2; };
        Atom_sp	getAtom3() { return this->_Atom3; };
	double	getT0()	{return this->term.t0;};
	double	getT();
	double	getTheta() {return this->getT();};

	double	getT0Degrees()	{return this->getT0()/0.0174533;};
	double	getTDegrees()	{return this->getT()/0.0174533;};
	double	getThetaDegrees(){return this->getTDegrees();};


	adapt::QDomNode_sp	asXml();
	void		parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml,
					AtomTable_sp	atomTable );
        void defineFrom( FFAngle_sp term, EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, double scale);
	void defineMissing( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3);
        core::List_sp encode() const;
        void decode(core::List_sp alist);
};
};

namespace translate {

template <>
struct	to_object<chem::EnergyAngle >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::EnergyAngle& angle)
  {
    return angle.encode();
  }
};

template <>
struct	from_object<chem::EnergyAngle>
{
  typedef	chem::EnergyAngle	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR("Implement me");
        }
};
};

namespace chem {

double	_evaluateEnergyOnly_Angle(
    int I1, int I2, int I3, core::T_sp activeAtomMask,
    num_real x1, num_real y1, num_real z1,
    num_real x2, num_real y2, num_real z2,
    num_real x3, num_real y3, num_real z3,
    num_real kt, num_real t0 );


FORWARD(EnergyAngle);

class EnergyAngle_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyAngle_O,"EnergyAngle",EnergyComponent_O);
public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
public: // virtual functions inherited from Object
    void	initialize();
public:
    typedef EnergyAngle	TermType;
public: // instance variables
    gctools::Vec0<TermType>	_Terms;
    gctools::Vec0<TermType>	_BeyondThresholdTerms;
public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
	
public:
    virtual size_t numberOfTerms() { return this->_Terms.size();};
    CL_DEFMETHOD core::T_mv safe_amber_energy_angle_term(size_t index) {
      if (index >= this->numberOfTerms() ) {
        SIMPLE_ERROR("Illegal term index {} must be less than {}" , index , this->_Terms.size() );
      }
      return Values(core::DoubleFloat_O::create(this->_Terms[index].term.kt),
                    core::DoubleFloat_O::create(this->_Terms[index].term.t0),
                    core::make_fixnum(this->_Terms[index].term.I1),
                    core::make_fixnum(this->_Terms[index].term.I2),
                    core::make_fixnum(this->_Terms[index].term.I3)
                    );
    }

public:	// Creation class functions



public:
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

  virtual	core::List_sp checkForBeyondThresholdInteractionsWithPosition(NVector_sp pos, double threshold );

    virtual string	beyondThresholdInteractionsAsString();

  void addAngleTerm(AtomTable_sp at, Atom_sp a1, Atom_sp a2, Atom_sp a3, double kt, double t0);


  core::List_sp lookupAngleTerms(AtomTable_sp at, Atom_sp a1, Atom_sp a2 , Atom_sp a3, core::HashTable_sp atomTypes );

  EnergyAngle_sp copyFilter(core::T_sp keepInteractionFactory);

    EnergyAngle_O( const EnergyAngle_O& ss ); //!< Copy constructor

  EnergyAngle_O() : EnergyComponent_O() {};
};





};

#endif //]
