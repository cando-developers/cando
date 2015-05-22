       
       

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
#include "core/common.h"
#include "addon/vector3.h"
#include "energyComponent.h"



#include "chemPackage.h"


namespace chem {
class EnergyAtom;

SMART(FFAngle);

inline	string	XmlTag_Angle() { return "Angle"; };
inline	string	XmlTag_EnergyAngle() { return "EnergyAngle"; };

typedef struct {
	REAL		kt;	//!< Angle bend force constant, this must match Mathematica code!
	REAL		t0;	//!< Angle bend equilibrium angle (Radians), this must match Mathematica code!
	INT		I1;	//!< i*3 index into coordinate vector for atom1, this must match Mathematica code!
	INT		I2;	//!< i*3 index into coordinate vector for atom2, this must match Mathematica code!
	INT		I3;	//!< i*3 index into coordinate vector for atom3, this must match Mathematica code!
} TermAngle;

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
#include "_Angle_debugEvalDeclares.cc"
#endif

public:
//	void	archive( core::ArchiveP node );

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

	addon::QDomNode_sp	asXml(core::Lisp_sp);
	void		parseFromXmlUsingAtomTable(addon::QDomNode_sp xml,
					AtomTable_sp	atomTable );
        void defineFrom( FFAngle_sp term, EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, double scale);
	void defineMissing( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3);
    EnergyAngle();
	virtual ~EnergyAngle();
};



double	_evaluateEnergyOnly_Angle(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3,
		double kt, double t0 );


class EnergyAngle_O : public EnergyComponent_O
{
    LISP_BASE1(EnergyComponent_O);
    LISP_CLASS(chem,ChemPkg,EnergyAngle_O,"EnergyAngle");
public:
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);

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
    virtual int numberOfTerms() { return this->_Terms.size();};

public:	// Creation class functions



public:
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




    EnergyAngle_O( const EnergyAngle_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyAngle_O);
};





};

TRANSLATE(chem::EnergyAngle_O);
#endif //]


