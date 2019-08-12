/*
    File: energyPointToLineRestraint.h
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
 *	energyPointToLineRestraint.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyPointToLineRestraint_H  //[
#define	EnergyPointToLineRestraint_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>



namespace chem
{


    struct TermPointToLineRestraint {
	double		ka;
	double		ra;
	int		I1;
	double		xa;
	double		ya;
	double		za;
	double		xb;
	double		yb;
	double		zb;
    };

/*! Store a PointToLineRestraint energy term
 */
class	EnergyPointToLineRestraint : public EnergyTerm {
public:
	string	className()	{ return "EnergyPointToLineRestraint"; };
public:
	TermPointToLineRestraint	term;
		// Variables
        Atom_sp      _Atom1;
		// Threshold
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_PointToLineRestraint_debugEvalDeclares.cc>
#endif
        Atom_sp	getAtom() { return this->_Atom1; };
	double		getXa() { return this->term.xa; };
	double		getYa() { return this->term.ya; };
	double		getZa() { return this->term.za; };
	double		getXb() { return this->term.xb; };
	double		getYb() { return this->term.yb; };
	double		getZb() { return this->term.zb; };
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );

    EnergyPointToLineRestraint();
	virtual ~EnergyPointToLineRestraint();

};



double	_evaluateEnergyOnly_PointToLineRestraint(
		double x1, double y1, double z1,
		double xa, double ya, double za,
		double xb, double yb, double zb,
		double ka, double ra );

class EnergyPointToLineRestraint_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyPointToLineRestraint_O,"EnergyPointToLineRestraint",EnergyComponent_O);
public:
public: // virtual functions inherited from Object
    void	initialize();
public:
    typedef EnergyPointToLineRestraint	TermType;
public: // instance variables
    gctools::Vec0<TermType>		_Terms;
    gctools::Vec0<TermType>		_BeyondThresholdTerms;

public:	// Creation class functions
public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
public:
    virtual size_t numberOfTerms() { return this->_Terms.size();};


public:
    void addTerm(const TermType& term);
    virtual void dumpTerms();

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

    // virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

    virtual string	beyondThresholdInteractionsAsString();


public:
    EnergyPointToLineRestraint_O( const EnergyPointToLineRestraint_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyPointToLineRestraint_O);
};




};

TRANSLATE(chem::EnergyPointToLineRestraint_O);
#endif //]
