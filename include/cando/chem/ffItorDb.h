/*
    File: ffItorDb.h
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
 *	ffItor.h
 *
 *	Maintains a database of itor terms
 */

#ifndef FFITOR_H
#define	FFITOR_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/symbolMap.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/ffBaseDb.h>
#include <clasp/core/smallMap.h>

namespace       chem
{

#define	FFItorWildCardCode	-1


    inline	string	XmlTag_FFItor() { return "FFItor";};

    SMART(FFItor);
    class FFItor_O : public FFParameter_O
    {
	LISP_CLASS(chem,ChemPkg,FFItor_O,"FFItor",FFParameter_O);
    public:
	static int const IMaxPeriodicity = 6;
    public:
	void initialize();
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
    public:
	core::Symbol_sp         _T1;
	core::Symbol_sp 	_T2;
	core::Symbol_sp 	_T3;
	core::Symbol_sp 	_T4;
	bool      	_hasPeriodicity[IMaxPeriodicity];
	double      _Vs_kJ[IMaxPeriodicity];
	double      _PhaseRads[IMaxPeriodicity];

	void	setTypes( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4);

	bool hasPeriodicity(int period) const;
	double  getV_kJ(int period) const;
	double  getV_kCal(int period) const;
	void    setV_kJ(int period, double v);
	void    setV_kCal(int period, double v);
	double  getPhaseRad(int period) const;
	void    setPhaseRad(int period, double v);

	void    mergeWith(FFItor_sp itor);


	virtual	string	levelDescription();
	virtual ParameterType type() { return itor; };
	DEFAULT_CTOR_DTOR(FFItor_O);
    };




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
    SMART(FFItorDb);
    class FFItorDb_O : public FFBaseDb_O
    {
	LISP_CLASS(chem,ChemPkg,FFItorDb_O,"FFItorDb",FFBaseDb_O);

    public:
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
    public:
	typedef	gctools::Vec0<FFItor_sp>::iterator	FFItor_spIterator;
        gctools::Vec0<FFItor_sp>		_Terms;
	adapt::SymbolMap<FFItor_O>		_Lookup;

	FFItor_spIterator	begin() { return this->_Terms.begin();};
	FFItor_spIterator	end() { return this->_Terms.end();};

	void	add( FFItor_sp itor );

        gc::Nilable<FFItor_sp> findExactTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );
	bool    hasExactTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );

	//! Look for exact term then general one
        gc::Nilable<FFItor_sp> findBestTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );
        bool    hasBestTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );

        void    cantFind(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 );

	//! Dump all ptors to stdout that match these types
	DEFAULT_CTOR_DTOR(FFItorDb_O);
    };


};

TRANSLATE(chem::FFItor_O);
TRANSLATE(chem::FFItorDb_O);
#endif
