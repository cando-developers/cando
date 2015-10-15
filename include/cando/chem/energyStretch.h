       
       

//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyStretch.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyStretch1_H  //[
#define	EnergyStretch1_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>

namespace       chem {

class EnergyAtom;
SMART(FFStretch);

inline	string	XmlTag_Stretch() { return "Stretch"; };
inline	string	XmlTag_EnergyStretch() { return "EnergyStretch"; };


struct TermStretch {
	REAL	kb;	//!< Stretch force constant, this must match Mathematica code!
	REAL	r0;	//!< Stretch equilibrium distance, this must match Mathematica code!
	INT	I1;	//!< i*3 index into coordinate vector for atom1, this must match Mathematica code!
	INT	I2;	//!< i*3 index into coordinate vector for atom2, this must match Mathematica code!
};


/*! Store a Stretch energy term
 */
class EnergyStretch : public EnergyTerm {
public:
	string	className()	{ return "EnergyStretch"; };
public:
	TermStretch	term;
        gc::Nilable<Atom_sp>		_Atom1;
        gc::Nilable<Atom_sp>		_Atom2;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_Stretch_debugEvalDeclares.cc>
#endif


	gc::Nilable<Atom_sp>	getAtom1() {return this->_Atom1;};
        gc::Nilable<Atom_sp>	getAtom2() {return this->_Atom2;};
	double	getR0()	{return this->term.r0;};
	double	getR();

public:
//	void	archive(core::ArchiveP node);
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );


	EnergyStretch();
	virtual ~EnergyStretch();
		//! Initialize self from a FFStretch parameter and two EnergyAtom(s)
	void defineFrom( FFStretch_sp stretch, EnergyAtom *ea1, EnergyAtom *ea2, double scale);
		//! Handle situation when FFStretch term could not be found for these two atoms
	void defineMissing( EnergyAtom *ea1, EnergyAtom *ea2);
//        void defineFrom( EnergyAtom *ea1, EnergyAtom *ea2);
};



double	_evaluateEnergyOnly_Stretch (
		double x1,
		double y1,
		double z1,
		double x2,
		double y2,
		double z2,
		double r0,
		double kb );






class EnergyStretch_O : public EnergyComponent_O
{
    LISP_BASE1(EnergyComponent_O);
    LISP_CLASS(chem,ChemPkg,EnergyStretch_O,"EnergyStretch");
public:
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

public:
    typedef EnergyStretch	TermType;
public: // instance variables
    gctools::Vec0<TermType>	_Terms;
    gctools::Vec0<TermType>	_BeyondThresholdTerms;
public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
	
public:
    virtual int numberOfTerms() { return this->_Terms.size();};
    void addTerm(const TermType& term);
    virtual void dumpTerms();

    virtual void setupHessianPreconditioner(NVector_sp nvPosition,
					    AbstractLargeSquareMatrix_sp m );
    virtual	void evaluateAll ( NVector_sp 	pos,
				   bool 		calcForce,
				   NVector_sp 	force,
				   bool		calcDiagonalHessian,
				   bool		calcOffDiagonalHessian,
				   AbstractLargeSquareMatrix_sp	hessian,
				   NVector_sp	hdvec,
				   NVector_sp dvec);

    virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
	NVector_sp pos );

    virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

    virtual string	beyondThresholdInteractionsAsString();


public:
    EnergyStretch_O( const EnergyStretch_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyStretch_O);
};







};

TRANSLATE(chem::EnergyStretch_O);
#endif //]


