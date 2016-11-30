/*
    File: forceField.h
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
 *	forceField.h
 *
 *	Maintains a parameter set for the Moe general force field
 *	from these parameters we can derive any force field parameters
 *
 */

#ifndef FORCEFIELD_H
#define	FORCEFIELD_H

#include <clasp/core/common.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
//#include	"conformation.h"
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/ffItorDb.h>
#include <cando/chem/ffPtorDb.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/ffVdwDb.h>



namespace       chem {


extern string	ForceField_XmlName;



typedef	enum	{ ffStr, ffAng, ffStb, ffOop, ffItor, ffTor, ffVdw, ffEle, ffSol } DisableEnum;


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
class	InfoDb_O;
    typedef	gctools::smart_ptr<InfoDb_O>	InfoDb_sp;
SMART(InfoDb);
class InfoDb_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,InfoDb_O,"InfoDb",core::CxxObject_O);

 public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
 public:
  adapt::SymbolMap<core::Str_O>	_database;
 public:
  void absorb(InfoDb_sp other) {
    for ( auto it : other->_database ) {
      this->addInfo(it.first,it.second);
    }
  };
  void	addInfo( core::Symbol_sp key, core::Str_sp data );
  DEFAULT_CTOR_DTOR(InfoDb_O);
};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

SMART(ForceField);
class ForceField_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,ForceField_O,"ForceField",core::CxxObject_O);

 public:
  void initialize();
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
 public:
//		/*! Read the forceField from a file*/
//	static ForceField_sp	open_ForceField(const string& fn);
  CL_LISPIFY_NAME("make-ForceField");
  CL_DEF_CLASS_METHOD static ForceField_sp make() { GC_ALLOCATE(ForceField_O,ff); return ff;};
 public:
  string				_Title;
  string				_Ref;
  gctools::Vec0<core::Symbol_sp>		_SingleBondMultiBondDistinctions;
  InfoDb_sp			_Info;
  FFTypesDb_sp	_Types;
  FFStretchDb_sp	_Stretches;
  FFAngleDb_sp	_Angles;
  FFItorDb_sp	_Itors;
  FFPtorDb_sp	_Ptors;
  FFNonbondDb_sp	_Nonbonds;
  FFVdwDb_sp		_Vdws;
//	bool				_Disable[sizeof(DisableEnum)];
#if 0
  RPFFEquivalenceDb		Equivalences;
  RPFFOutOfPlaneDb		OutOfPlanes;
  FFStretch_spBendDb		StretchBends;
#endif
 public:
 ForceField_O() :
  _Info(_Unbound<InfoDb_O>()),
    _Types(_Unbound<FFTypesDb_O>()),
    _Stretches(_Unbound<FFStretchDb_O>()),
    _Angles(_Unbound<FFAngleDb_O>()),
    _Itors(_Unbound<FFItorDb_O>()),
    _Ptors(_Unbound<FFPtorDb_O>()),
    _Nonbonds(_Unbound<FFNonbondDb_O>()),
    _Vdws(_Unbound<FFVdwDb_O>())
    {};
          
  void	addSingleBondMultiBondDistinctionType(core::Symbol_sp s) {
    this->_SingleBondMultiBondDistinctions.push_back(s);
  };

  CL_NAME("getTypes");
  CL_DEFMETHOD 	FFTypesDb_sp	getTypes() { return this->_Types; };

  CL_NAME("getStretchDb");
  CL_DEFMETHOD 	FFStretchDb_sp getStretchDb() { return this->_Stretches;};
  CL_NAME("getAngleDb");
  CL_DEFMETHOD 	FFAngleDb_sp getAngleDb() { return this->_Angles;};
  CL_NAME("getItorDb");
  CL_DEFMETHOD 	FFItorDb_sp getItorDb() { return this->_Itors;};
  CL_NAME("getPtorDb");
  CL_DEFMETHOD 	FFPtorDb_sp getPtorDb() { return this->_Ptors;};
  CL_NAME("getNonbondDb");
  CL_DEFMETHOD 	FFNonbondDb_sp	getNonbondDb() { return this->_Nonbonds; };
  CL_NAME("getVdwDb");
  CL_DEFMETHOD 	FFVdwDb_sp getVdwDb() { return this->_Vdws;};

  void absorb(ForceField_sp other);
  void	assignTypes(Matter_sp matter);
  void	setTitle(const string& title);
  void	setInfoDb( InfoDb_sp Info );
  void	setFFTypeDb( FFTypesDb_sp Types);
  void	setFFStretchDb( FFStretchDb_sp Stretches);
  void	setFFAngleDb( FFAngleDb_sp Angles);
  void	setFFItorDb( FFItorDb_sp Itors);
  void	setFFPtorDb( FFPtorDb_sp Ptors);
  void	setFFNonbondDb(FFNonbondDb_sp Nonbonds );
  void	setFFVdwDb(FFVdwDb_sp Vdws );

//	void		parseFromMoeStream(istream in);
  void		parseFromMoeFormatFileName(string name);


  void	pointAllPartsToForceField(ForceField_sp ff);

  void	saveAs(const string& fileName);
};



//
// Some forcefield values may be undefined
// they will be set to UNDEFINED_NUMBER
// test this with isUndefined(x)
//
#define	UNDEFINED_NUMBER	9.0e99
#define	UNDEFINED_CUTOFF	9.0e98
inline	bool	isUndefined(double x) {return (x)>=UNDEFINED_CUTOFF;}



};


TRANSLATE(chem::InfoDb_O);
TRANSLATE(chem::ForceField_O);
#endif
