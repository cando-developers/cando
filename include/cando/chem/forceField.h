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
    gc::SmallMap<core::Symbol_sp,core::String_sp>	_database;
  public:
    void forceFieldMerge(InfoDb_sp other) {
      for ( auto it : other->_database ) {
        this->addInfo(it.first,it.second);
      }
    };
    void	addInfo( core::Symbol_sp key, core::String_sp data );
    DEFAULT_CTOR_DTOR(InfoDb_O);
  };

};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

 

template <>
struct gctools::GCInfo<chem::ForceField_O> {
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {
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
  CL_LISPIFY_NAME("ForceField/make");
  CL_DEF_CLASS_METHOD static ForceField_sp make() {
    auto ff = gctools::GC<ForceField_O>::allocate_with_default_constructor();
    return ff;
  };
public:
  core::SimpleBaseString_sp      _title;
  gctools::Vec0<core::Symbol_sp> _SingleBondMultiBondDistinctions;
  InfoDb_sp			 _Info;
  FFTypesDb_sp	                 _Types;
  FFStretchDb_sp	         _Stretches;
  FFAngleDb_sp	                 _Angles;
  FFItorDb_sp	                 _Itors;
  FFPtorDb_sp	                 _Ptors;
  FFNonbondDb_sp	         _Nonbonds;
  FFLKSolvationDb_sp	         _LKSolvation;
  FFVdwDb_sp		         _Vdws;
public:
  ForceField_O() :
    _Info(unbound<InfoDb_O>()),
    _Types(unbound<FFTypesDb_O>()),
    _Stretches(unbound<FFStretchDb_O>()),
    _Angles(unbound<FFAngleDb_O>()),
    _Itors(unbound<FFItorDb_O>()),
    _Ptors(unbound<FFPtorDb_O>()),
    _Nonbonds(unbound<FFNonbondDb_O>()),
    _LKSolvation(unbound<FFLKSolvationDb_O>()),
    _Vdws(unbound<FFVdwDb_O>())
  {};
          
  void	addSingleBondMultiBondDistinctionType(core::Symbol_sp s) {
    this->_SingleBondMultiBondDistinctions.push_back(s);
  };

  core::SimpleBaseString_sp getTitle() const;
  
  CL_LISPIFY_NAME("getTypes");
  CL_DEFMETHOD FFTypesDb_sp getTypes() { return this->_Types; };

  CL_LISPIFY_NAME("getStretchDb");
  CL_DEFMETHOD 	FFStretchDb_sp getStretchDb() { return this->_Stretches;};
  CL_LISPIFY_NAME("getAngleDb");
  CL_DEFMETHOD 	FFAngleDb_sp getAngleDb() { return this->_Angles;};
  CL_LISPIFY_NAME("getItorDb");
  CL_DEFMETHOD 	FFItorDb_sp getItorDb() { return this->_Itors;};
  CL_LISPIFY_NAME("getPtorDb");
  CL_DEFMETHOD 	FFPtorDb_sp getPtorDb() { return this->_Ptors;};
  CL_LISPIFY_NAME("getNonbondDb");
  CL_DEFMETHOD 	FFNonbondDb_sp	getNonbondDb() { return this->_Nonbonds; };
  CL_LISPIFY_NAME("getLKSolvationDb");
  CL_DEFMETHOD 	FFLKSolvationDb_sp	getLKSolvationDb() { return this->_LKSolvation; };
  CL_LISPIFY_NAME("getVdwDb");
  CL_DEFMETHOD 	FFVdwDb_sp getVdwDb() { return this->_Vdws;};

  void forceFieldMerge(ForceField_sp other);
  void	assignTypes(Matter_sp matter, core::HashTable_sp atom_types);
  void	setTitle(const string& title);
  void	setInfoDb( InfoDb_sp Info );
  void	setFFTypeDb( FFTypesDb_sp Types);
   void	setFFStretchDb( FFStretchDb_sp Stretches);
   void	setFFAngleDb( FFAngleDb_sp Angles);
   void	setFFItorDb( FFItorDb_sp Itors);
   void	setFFPtorDb( FFPtorDb_sp Ptors);
   void	setFFNonbondDb(FFNonbondDb_sp Nonbonds );
   void	setFFLKSolvationDb(FFLKSolvationDb_sp Nonbonds );
   void	setFFVdwDb(FFVdwDb_sp Vdws );

//	void		parseFromMoeStream(istream in);
  void		parseFromMoeFormatFileName(string name);

  void	pointAllPartsToForceField(ForceField_sp ff);

  string __repr__() const;

  //  void writeAmberForceFieldParameterFile(core::T_sp stream, core::T_sp system);
};



//
// Some forcefield values may be undefined
// they will be set to UNDEFINED_NUMBER
// test this with isUndefined(x)
//
#define	UNDEFINED_NUMBER	9.0e99
#define	UNDEFINED_CUTOFF	9.0e98
inline	bool	isUndefined(double x) {return (x)>=UNDEFINED_CUTOFF;}


SMART(CombinedForceField);
class CombinedForceField_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,CombinedForceField_O,"CombinedForceField",core::CxxObject_O);
public:
  core::List_sp _ForceFields;
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
  static CombinedForceField_sp make();
public:
  void addShadowingForceField(ForceField_sp forceField, core::T_sp info);
  core::List_sp forceFieldsAsList() const;
  void clear();

  core::T_sp assignForceFieldTypes(Matter_sp molecule, core::HashTable_sp atom_types );
  void assignMolecularEnergyTables(Matter_sp molecule,core::T_sp keepInteractionFactory);
  
  CombinedForceField_O() : _ForceFields(nil<core::T_O>()) {};
};

};

#endif


