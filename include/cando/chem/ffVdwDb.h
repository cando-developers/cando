/*
    File: ffVdwDb.h
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
