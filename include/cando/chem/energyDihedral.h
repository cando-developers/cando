/*
    File: energyDihedral.h
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
 *	energyDihedral.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyDihedral_H  //[
#define	EnergyDihedral_H
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

SMART(FFPtor);
SMART(FFItor);

inline	string	XmlTag_Dihedral() { return "Dihedral"; };
inline	string	XmlTag_EnergyDihedral() { return "EnergyDihedral"; };

struct TermDihedral
{
	REAL		sinPhase;	//!< Sine of phase, this must match Mathematica code!
	REAL		cosPhase;	//!< Cosine of phase, this must match Mathematica code!
	REAL		V;	//!< Dihedral force constant, this must match Mathematica code!
        REAL          	DN;	//!< Dihedral multiplicity as double, this must match Mathematica code!
        INT		IN;	//!< Dihedral multiplicity as int, this must match Mathematica code!
	INT		I1;	//!< i*3 index into coordinate array, this must match Mathematica code!
	INT		I2;	//!< i*3 index into coordinate array, this must match Mathematica code!
	INT		I3;	//!< i*3 index into coordinate array, this must match Mathematica code!
	INT		I4;	//!< i*3 index into coordinate array, this must match Mathematica code!
  TermDihedral(REAL sinp, REAL cosp, REAL v, REAL dn, INT in, INT i1, INT i2, INT i3, INT i4)
    : sinPhase(sinp),
      cosPhase(cosp),
      V(v),
      DN(dn),
      IN(in),
      I1(i1), I2(i2), I3(i3), I4(i4) {};
  TermDihedral() {};
};

/*! Dihedral energy term
 */
class EnergyDihedral : public EnergyTerm 
{
public:
	string	className()	{ return "EnergyDihedral"; };
public:
                // Parameters
//	string		_Type1;
//	string		_Type2;
//	string		_Type3;
//	string		_Type4;
	bool		_Proper;
                // Variables
        Atom_sp      _Atom1;
        Atom_sp      _Atom2;
        Atom_sp      _Atom3;
        Atom_sp      _Atom4;
	REAL		_PhaseRad;
	TermDihedral	term;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_Dihedral_debugEvalDeclares.cc>
#endif
	REAL		_CalculatedDihedralDeviation;

        Atom_sp		getAtom1() { return this->_Atom1; };
        Atom_sp		getAtom2() { return this->_Atom2; };
        Atom_sp		getAtom3() { return this->_Atom3; };
        Atom_sp		getAtom4() { return this->_Atom4; };
	double		getDihedral();
	double		getDihedralDeviation();
	double		getCalculatedDihedralDeviation() { return this->_CalculatedDihedralDeviation;};

public:
 EnergyDihedral(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double phase, bool proper, size_t i1, size_t i2, size_t i3, size_t i4, double sinp, double cosp, double v, double dn, int in) : term(sinp,cosp,v,dn,in,i1,i2,i3,i4),_Atom1(a1),_Atom2(a2),_Atom3(a3),_Atom4(a4), _PhaseRad(phase), _Proper(proper) {};
  EnergyDihedral() {};
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );

        void defineFrom( int n, FFPtor_sp term, EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale);
	void defineMissingProper( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4);
        void defineFrom( int n, FFItor_sp term, EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale);
        core::List_sp encode() const;
        void decode(core::List_sp alist);
};
};



namespace translate {

template <>
struct	to_object<chem::EnergyDihedral >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::EnergyDihedral& dihedral)
  {
    return dihedral.encode();
  }
};

template <>
struct	from_object<chem::EnergyDihedral>
{
  typedef	chem::EnergyDihedral	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR("Implement me");
        }
};
};

namespace chem {

double  _evaluateEnergyOnly_Dihedral(
                double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3,
                double x4, double y4, double z4,
                double V, double DN, int IN,
                double cosPhase, double sinPhase );



FORWARD(EnergyDihedral);

class EnergyDihedral_O : public EnergyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyDihedral_O,"EnergyDihedral",EnergyComponent_O);
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public: // virtual functions inherited from Object
  void	initialize();

public:
  typedef EnergyDihedral	TermType;
public: // instance variables
  gctools::Vec0<TermType>	_Terms;
  gctools::Vec0<TermType>	_BeyondThresholdTerms;

public:	// Creation class functions
  typedef gctools::Vec0<TermType>::iterator iterator;
  iterator begin() { return this->_Terms.begin(); };
  iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
public:
  virtual size_t numberOfTerms() { return this->_Terms.size();};
    

public:
  void addTerm(const TermType& term);
  virtual void dumpTerms(core::HashTable_sp atomTypes);

  CL_DEFMETHOD core::T_mv safe_amber_energy_dihedral_term(size_t index) {
    if (index >= this->numberOfTerms() ) {
      SIMPLE_ERROR("Illegal term index {} must be less than {}" , index , this->_Terms.size() );
    }
    return Values(core::DoubleFloat_O::create(this->_Terms[index].term.cosPhase), // Note switched order
                  core::DoubleFloat_O::create(this->_Terms[index].term.sinPhase),
                  core::DoubleFloat_O::create(this->_Terms[index].term.V),
                  core::DoubleFloat_O::create(this->_Terms[index].term.DN),
                  core::make_fixnum(this->_Terms[index].term.I1),
                  core::make_fixnum(this->_Terms[index].term.I2),
                  core::make_fixnum(this->_Terms[index].term.I3),
                  core::make_fixnum(this->_Terms[index].term.I4)
                  );
  }

  virtual core::List_sp extract_vectors_as_alist() const;
  virtual void fill_from_vectors_in_alist(core::List_sp vectors);

    
  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m );
  virtual double evaluateAllComponent( ScoringFunction_sp scorer,
                                       NVector_sp 	pos,
                                       core::T_sp componentEnergy,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec,
                                       gc::Nilable<NVector_sp> dvec);

  virtual double evaluateAllComponentSingle(
      gctools::Vec0<EnergyDihedral>::iterator di_start,
      gctools::Vec0<EnergyDihedral>::iterator di_end,
      ScoringFunction_sp scorer,
      NVector_sp 	pos,
      bool 		calcForce,
      gc::Nilable<NVector_sp> 	force,
      bool		calcDiagonalHessian,
      bool		calcOffDiagonalHessian,
      gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
      gc::Nilable<NVector_sp>	hdvec,
      gc::Nilable<NVector_sp> dvec);

  virtual double evaluateAllComponentSimd8(
      gctools::Vec0<EnergyDihedral>::iterator di_start8,
      gctools::Vec0<EnergyDihedral>::iterator di_end8,
      ScoringFunction_sp scorer,
      NVector_sp 	pos,
      bool 		calcForce,
      gc::Nilable<NVector_sp> 	force,
      bool		calcDiagonalHessian,
      bool		calcOffDiagonalHessian,
      gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
      gc::Nilable<NVector_sp>	hdvec,
      gc::Nilable<NVector_sp> dvec);

  virtual double evaluateAllComponentSimd4(
      gctools::Vec0<EnergyDihedral>::iterator di_start4,
      gctools::Vec0<EnergyDihedral>::iterator di_end4,
      ScoringFunction_sp scorer,
      NVector_sp 	pos,
      bool 		calcForce,
      gc::Nilable<NVector_sp> 	force,
      bool		calcDiagonalHessian,
      bool		calcOffDiagonalHessian,
      gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
      gc::Nilable<NVector_sp>	hdvec,
      gc::Nilable<NVector_sp> dvec);

  virtual double evaluateAllComponentSimd2(
      gctools::Vec0<EnergyDihedral>::iterator di_start2,
      gctools::Vec0<EnergyDihedral>::iterator di_end2,
      ScoringFunction_sp scorer,
      NVector_sp 	pos,
      bool 		calcForce,
      gc::Nilable<NVector_sp> 	force,
      bool		calcDiagonalHessian,
      bool		calcOffDiagonalHessian,
      gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
      gc::Nilable<NVector_sp>	hdvec,
      gc::Nilable<NVector_sp> dvec);

  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
      NVector_sp pos );

    // virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

  virtual string	beyondThresholdInteractionsAsString();


  void addDihedralTerm(AtomTable_sp at, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double phase, bool proper, double v, int multiplicity);

  core::List_sp lookupDihedralTerms(AtomTable_sp at, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, core::HashTable_sp atomTypes );

  EnergyDihedral_sp copyFilter(core::T_sp keepInteraction);

public:
  EnergyDihedral_O( const EnergyDihedral_O& ss ); //!< Copy constructor

  EnergyDihedral_O() {};
};


};

TRANSLATE(chem::EnergyDihedral_O);
#endif //]
