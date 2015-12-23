       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffVdw.h
 *
 *	Maintains a database of vdw terms
 */

#ifndef FFVDW_H
#define	FFVDW_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
//#include	"conformation.h"
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/ffBaseDb.h>


#include <cando/chem/chemPackage.h>


namespace chem {



//#define	FFVdwWildCardCode	-1

class	FFVdwDb_O;
    typedef	gctools::smart_ptr<FFVdwDb_O>	FFVdwDb_sp;

class	FFVdw_O;
    typedef	gctools::smart_ptr<FFVdw_O>	FFVdw_sp;
SMART(FFVdw);
class FFVdw_O : public FFParameter_O
{
    LISP_CLASS(chem,ChemPkg,FFVdw_O,"FFVdw",FFParameter_O);

public:
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
public:
	string		_T1;
	string		_T2;
	double		_R;
	double		_Eps;
	double		_M;
	double		_N;
	double		_Edep;
//virtual	string	levelDescription();
virtual ParameterType type() { return vdw; };
DEFAULT_CTOR_DTOR(FFVdw_O);
};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
SMART(FFVdwDb);
class FFVdwDb_O : public FFBaseDb_O
{
    LISP_CLASS(chem,ChemPkg,FFVdwDb_O,"FFVdwDb",FFBaseDb_O);

public:
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
public:
    gctools::Vec0<FFVdw_sp>	_Terms;

    void	add( FFVdw_sp vdw ) {SUBIMP();};

	DEFAULT_CTOR_DTOR(FFVdwDb_O);
};



};


TRANSLATE(chem::FFVdw_O);
TRANSLATE(chem::FFVdwDb_O);
#endif
