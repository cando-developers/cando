/*
    File: ffStretchDb.h
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
 *	ffStretch.h
 *
 *	Maintains a database of stretch types
 */

#ifndef FFSTR_H
#define	FFSTR_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/adapt/quickDom.h>
#include <cando/adapt/symbolMap.h>
#include <clasp/core/smallMap.h>
#include <cando/chem/ffBaseDb.h>



namespace       chem {



#define	FFStretchWildCardCode	-1

    SMART(FFStretch);
    class FFStretch_O : public FFParameter_O
    {
	LISP_CLASS(chem,ChemPkg,FFStretch_O,"FFStretch",FFParameter_O);
    public:
	void initialize();
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
    public:
        core::Symbol_sp	_Type1;
        core::Symbol_sp	_Type2;
	double		_R0_Nanometer;
	double		_Kb_kjPerNanometerSquared;
//	double		_K3;
//	double		_K4;
//	double		_Bci;

    public:
        static FFStretch_sp create_missing(core::Symbol_sp t1, core::Symbol_sp t2) {
          FFStretch_sp stretch = FFStretch_O::create();
          stretch->_Level = unknown;
          stretch->_Type1 = t1;
          stretch->_Type2 = t2;
          return stretch;
        }
    public:
	void setKb_kjPerNanometerSquared(double val);
	void setR0_Nanometer(double val);

	double getR0_Nanometer() const;
	double getR0_Angstrom() const;
	double getKb_kjPerNanometerSquared() const;
	double getKb_kcalPerAngstromSquared() const;
	virtual	string	levelDescription();
	DEFAULT_CTOR_DTOR(FFStretch_O);
    };




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
    class	EstimateStretch {
    public:
      core::Symbol_sp		_ti;
      core::Symbol_sp		_tj;
      double		_rij;
      double		_lnKij;
    };


    class	FFStretchDb_O;
    SMART(FFStretchDb);
    class FFStretchDb_O : public FFParameterBaseDb_O
    {
      LISP_CLASS(chem,ChemPkg,FFStretchDb_O,"FFStretchDb",FFParameterBaseDb_O);
    public:
      static FFStretchDb_sp make();
      FFStretchDb_O() {};
    public:
      bool fieldsp() const { return true; };
      void	fields(core::Record_sp node);
    public:
      gctools::SmallMap<core::Symbol_sp,EstimateStretch>	_EstimateStretch; // Not archived

      void	clearEstimateStretch();
      void	addEstimateStretch(core::Symbol_sp ti, core::Symbol_sp tj, double rij, double lnKij );
      void	_addEstimateStretch(const EstimateStretch& es);
      void	add( FFStretch_sp str );
      core::T_sp findTermForTypes(core::Symbol_sp a1, core::Symbol_sp a2 );

      void forceFieldMerge(FFBaseDb_sp other);
    };



};



TRANSLATE(chem::FFStretch_O);
TRANSLATE(chem::FFStretchDb_O);
#endif
