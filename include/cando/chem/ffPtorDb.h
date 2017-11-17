/*
    File: ffPtorDb.h
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
 *	ffPtor.h
 *
 *	Maintains a database of ptor terms
 */

#ifndef FFPTOR_H
#define	FFPTOR_H
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


#include <cando/chem/chemPackage.h>


namespace chem
{


#define	FFPtorWildCardCode	-1


    SMART(FFPtor);
    class FFPtor_O : public FFParameter_O
    {
	LISP_CLASS(chem,ChemPkg,FFPtor_O,"FFPtor",FFParameter_O);
    public:
	static int const MaxPeriodicity = 8;
    public:
	void initialize();
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
    public:
	core::Symbol_sp	_T1;
	core::Symbol_sp _T2;
	core::Symbol_sp _T3;
	core::Symbol_sp _T4;
    public:
        static FFPtor_sp create_missing(core::Symbol_sp t1,
                                 core::Symbol_sp t2,
                                 core::Symbol_sp t3,
                                 core::Symbol_sp t4)
        {
          FFPtor_sp ptor = FFPtor_O::create();
          ptor->_Level = unknown;
          ptor->_T1 = t1;
          ptor->_T2 = t2;
          ptor->_T3 = t3;
          ptor->_T4 = t4;
          return ptor;
        }
          
    public:
	bool		_HasPeriodicity[MaxPeriodicity];	// 6 terms can be stored
	double		_Vs_kJ[MaxPeriodicity];		// 6 V terms
        double          _PhaseRads[MaxPeriodicity]; 	// 6 phase terms
        bool            _HasScee;
        double          _Scee;
        bool            _HasScnb;
        double          _Scnb;
        string          _Comment;

	void	setTypes( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4);
	int	maxPeriodicity() const { return MaxPeriodicity; };
        bool	hasPeriodicity(int period) const;
        double  getV_kJ(int period) const;
        void    setV_kJ(int period, double v);
	double	getV_kCal(int period) const;
	void	setV_kCal(int period, double v);

        double  getPhaseRad(int period) const;
        void    setPhaseRad(int period, double pr);

        CL_DEFMETHOD void    setScee(double s) { this->_HasScee = true; this->_Scee = s;}
        CL_DEFMETHOD void    setScnb(double s) { this->_HasScnb = true; this->_Scnb = s;}
        CL_DEFMETHOD void    setHasScee(bool b) { this->_HasScee = b; }; 
        CL_DEFMETHOD void    setHasScnb(bool b) { this->_HasScnb = b; };
        CL_DEFMETHOD double  getScee() const { return this->_Scee; };
        CL_DEFMETHOD double  getScnb() const { return this->_Scnb; };
        CL_DEFMETHOD bool    getHasScee() const { return this->_HasScee;};
        CL_DEFMETHOD bool    getHasScnb() const { return this->_HasScnb;};

        CL_DEFMETHOD string  ptorComment() const { return this->_Comment; };
        CL_DEFMETHOD void    setPtorComment(const string& s) { this->_Comment = s; };
        
       void    mergeWith(FFPtor_sp ptor);

	virtual	ParameterType	type() { return ptor;};
	virtual	string		levelDescription();
    FFPtor_O() : _HasScee(false), _Scee(-1.0), _HasScnb(false), _Scnb(-1.0), _Comment("") {};
    };




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
    SMART(FFPtorDb);
    class FFPtorDb_O : public FFParameterBaseDb_O
    {
	LISP_CLASS(chem,ChemPkg,FFPtorDb_O,"FFPtorDb",FFParameterBaseDb_O);

    public:
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
    public:
	void	setTermFormat( core::VectorStrings vterms );
	void	add( FFPtor_sp ptor );
	void	addFFPtor( FFPtor_sp ptor );
        core::T_sp findExactTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );

	//! Look for exact term then general one
        core::T_sp findBestTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );
        void    cantFind(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 );
	void		initialize();
    void forceFieldMerge(FFBaseDb_sp bother);
	DEFAULT_CTOR_DTOR(FFPtorDb_O);
    };


};
#endif
