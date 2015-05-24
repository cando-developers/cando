       
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
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
//#include "geom/render.fwd.h"// energyNonbond.h wants DisplayList needs render.fwd.h

namespace       chem {

SMART(EnergyNonbond);
SMART(ForceField);
//SMART(DisplayList);
SMART(AtomTable);

typedef	struct
{
	double		dQ1Q2;
	double		dA;
	double		dC;
	int		I1; //!< i*3 index into coordinate vector, must match Mathematica code!
	int		I2; //!< i*3 index into coordinate vector, must match Mathematica code!
} TermNonBond;


inline	string	XmlTag_Nonbond() { return "Nonbond"; };
inline	string	XmlTag_EnergyNonbond() { return "EnergyNonbond"; };

/*! A Nonbond term, this is subclassed in AmberFunction to deal with Amber specific stuff.
 * I should combine energyFunction.cc and amberFunction.cc together.
 */
class EnergyNonbond : public EnergyTerm {
public:
	string	className()	{ return "EnergyNonbond"; };
public:
                // Parameters
                // Variables
	bool		_Is14;
        Atom_sp      	_Atom1;
        Atom_sp      	_Atom2;
	double		_A;
	double		_C;
	double		_Charge1;
	double		_Charge2;
	double		_RStar;
	TermNonBond	term;
//	double		_NonBondScale; // Used in Amber
//	double		_ElectrostaticScale; // Used in Amber
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/_Nonbond_debugEvalDeclares.cc>
#endif

	Atom_sp	getAtom1() { return this->_Atom1; };
	Atom_sp	getAtom2() { return this->_Atom2; };
	double	getDistance();
	bool	defineFrom(ForceField_sp	forceField,
					bool		is14,
    					EnergyAtom	*iea1,
    					EnergyAtom	*iea2,
					EnergyNonbond_sp nb);

public:
//	void	archive(core::ArchiveP node);
public:
	geom::QDomNode_sp	asXml(core::Lisp_sp);
	void	parseFromXmlRelativeToContainer( geom::QDomNode_sp xml, Matter_sp parent );
	void	parseFromXmlUsingAtomTable(geom::QDomNode_sp xml, AtomTable_sp atomTable );

public:
    EnergyNonbond();

};



double	_evaluateEnergyOnly_Nonbond(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double dA, double dC, double dQ1Q2 );


class EnergyNonbond_O : public EnergyComponent_O
{
    LISP_BASE1(EnergyComponent_O);
    LISP_CLASS(chem,ChemPkg,EnergyNonbond_O,"EnergyNonbond");
public:
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

public:
    typedef EnergyNonbond TermType;
public: // instance variables
    double		_DielectricConstant;
    double		_ScaleVdw;
    double		_ScaleElectrostatic;
    double		_EnergyVdw;
    double		_EnergyElectrostatic;
    gctools::Vec0<TermType>	_Terms;
    gctools::Vec0<TermType>	_BeyondThresholdTerms;

public:	// Creation class functions

public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
public:
    virtual int numberOfTerms() { return this->_Terms.size();};


public:
    void	setDielectricConstant(double d) { this->_DielectricConstant = d; };
    double	getDielectricConstant() { return this->_DielectricConstant; };

    void	setVdwScale(double d) { this->_ScaleVdw = d; };
    double	getVdwScale()	{return this->_ScaleVdw; };
    void	setElectrostaticScale(double d) { this->_ScaleElectrostatic = d; };
    double	getElectrostaticScale()	{return this->_ScaleElectrostatic; };

    double	getVdwEnergy() { return this->_EnergyVdw; };
    double	getElectrostaticEnergy() { return this->_EnergyElectrostatic; };


public:

    void addTerm(const TermType& term);
    virtual void dumpTerms();

    virtual	void	zeroEnergy();
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

//    int countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn, core::Lisp_sp );

    virtual	double	getEnergy();


    void constructFromAtomTable(AtomTable_sp atomTable, ForceField_sp forceField);

    

public:
    EnergyNonbond_O( const EnergyNonbond_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyNonbond_O);
};

};

TRANSLATE(chem::EnergyNonbond_O);
#endif //]


