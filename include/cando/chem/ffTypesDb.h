       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffTypes.h
 *
 *	Maintains databases and structures to store types and type assignement
 *	rules.
 */

#ifndef FFTYPES_H
#define	FFTYPES_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/chemInfo.h>
#include <cando/chem/ffBaseDb.h>


#include <cando/chem/chemPackage.h>


namespace chem {

inline	string	XmlTag_FFTypesDb()	{ return "FFTypesDb";};


///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
class	FFTypesDb_O;
    typedef	gctools::smart_ptr<FFTypesDb_O>	FFTypesDb_sp;
SMART(FFTypesDb);
class FFTypesDb_O : public FFBaseDb_O
{
    LISP_BASE1(FFBaseDb_O);
    LISP_CLASS(chem,ChemPkg,FFTypesDb_O,"FFTypesDb");

public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
public:
    gctools::Vec0<ChemInfo_sp>	_TypeAssignmentRules;

    void	add( chem::ChemInfo_sp ci ) {
	this->_TypeAssignmentRules.push_back(ci);
    }

    int	numberOfRules() { return this->_TypeAssignmentRules.size();};
    chem::ChemInfo_sp	getRule(uint index);

    void    assignTypes( chem::Matter_sp matter );
    void	initialize();


    DEFAULT_CTOR_DTOR(FFTypesDb_O);
};




};


TRANSLATE(chem::FFTypesDb_O);
#endif
