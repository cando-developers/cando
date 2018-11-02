/*
    File: ffAngleDb.h
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
 *	ffAngle.h
 *
 *	Maintains a database of angle terms
 */

#ifndef FFANGLE_H
#define	FFANGLE_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/symbolMap.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/ffBaseDb.h>
#include <clasp/core/smallMap.h>
#include <cando/chem/chemPackage.h>


namespace chem {

  FORWARD(FFStretchDb);
  
#define		FFAngleWildCardCode	-1

inline	string	XmlTag_FFAngle() { return "FFAngle";};

class	FFAngle_O;
    typedef	gctools::smart_ptr<FFAngle_O>	FFAngle_sp;
SMART(FFAngle);
class FFAngle_O : public FFParameter_O
{
    LISP_CLASS(chem,ChemPkg,FFAngle_O,"FFAngle",FFParameter_O);

public:
	void initialize();
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
 public:
        static FFAngle_sp create_missing(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3) {
          FFAngle_sp angle = FFAngle_O::create();
          angle->_Type1 = t1;
          angle->_Type2 = t2;
          angle->_Type3 = t3;
          angle->_Level = unknown;
          return angle;
        }

public:
    core::Symbol_sp 		_Type1;
    core::Symbol_sp 		_Type2;
    core::Symbol_sp 		_Type3;
	double		_AngRad;
	double		_K2__kJPerRadianSquared;
//	double		_K3;
//	double		_K4;
//	double		_Ub_k;
//	double		_Ub_len;

	double getAngle_Radian() const { return this->_AngRad;};

	void setK2_kCalPerRadianSquared(double val);
	void setK2_kJPerRadianSquared(double val);
	
	double getK2_kCalPerRadianSquared() const;
	double getK2_kJPerRadianSquared() const;

virtual	string	levelDescription();
virtual ParameterType type() { return angle; };

DEFAULT_CTOR_DTOR(FFAngle_O);
};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


inline	string	XmlTag_FFAngleDb() { return "FFAngleDb";};
class	FFAngleDb_O;
    typedef	gctools::smart_ptr<FFAngleDb_O>	FFAngleDb_sp;
SMART(FFAngleDb);
class FFAngleDb_O : public FFParameterBaseDb_O
{
    LISP_CLASS(chem,ChemPkg,FFAngleDb_O,"FFAngleDb",FFParameterBaseDb_O);

public:
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
public:
	string			_AngleFunction;
  gc::SmallMap<core::Symbol_sp,core::DoubleFloat_sp>		_ZConstants;	//<! Angle parameter constants from Wang et al. J. Comput. Chem 25, 1157-1174 (2004)
  gc::SmallMap<core::Symbol_sp, core::DoubleFloat_sp>		_CConstants;	//<! Angle parameter constants from Wang et al. J. Comput. Chem 25, 1157-1174 (2004)

public:
	void	setAngleFunction(const string& st) {this->_AngleFunction = st;};
	void	add( FFAngle_sp ang );

    void	addZConstant(core::Symbol_sp element, double value ) { this->_ZConstants.set(element,core::DoubleFloat_O::create(value));};
    void	addCConstant(core::Symbol_sp element, double value ) { this->_CConstants.set(element,core::DoubleFloat_O::create(value));};

    FFAngle_sp	findTerm( FFStretchDb_sp ffstretch, chem::Atom_sp a1, chem::Atom_sp a2, chem::Atom_sp a3);
	FFAngle_sp	estimateTerm( FFStretchDb_sp ffstretch, chem::Atom_sp a1, chem::Atom_sp a2, chem::Atom_sp a3 );

	void	initialize();

        void forceFieldMerge(FFBaseDb_sp bother);

	DEFAULT_CTOR_DTOR(FFAngleDb_O);
};

};

TRANSLATE(chem::FFAngle_O);
TRANSLATE(chem::FFAngleDb_O);
#endif
