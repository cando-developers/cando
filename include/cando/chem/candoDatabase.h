/*
    File: candoDatabase.h
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
 	candoDatabase.h
 */

#ifndef	candoDatabase_H
#define	candoDatabase_H


#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <clasp/core/common.h>
#include <cando/adapt/stringList.fwd.h>
#include <cando/geom/vector3.h>
#include <cando/chem/chemPackage.h>


namespace chem {



SMART(Aggregate);
SMART(Entity);
SMART(SpecificContextSet);
SMART(ScorerVirtualMachine);
FORWARD(Residue);
SMART(RepresentedEntityNameSet);
SMART(MonomerPack);
SMART(MonomerContext);
SMART(FrameRecognizer);
SMART(Constitution);
SMART(EntityNameSet);
FORWARD(MonomerCoordinates);

extern	string	ResourceConstitution;
	/*! Used to lookup entries in the Group, Constitution, Stereoisomer name space */
extern	string	ResourceEntity;
extern	string	ResourceFrameRecognizer;


SMART(CandoDatabase);

/*! A CandoDatabase stores Constitutions
	as a Monomer_spTypeDatabase.
	The Constitutions define everything you need
	to define the topology of any molecule.

*/
SMART(CandoDatabase);
class CandoDatabase_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,CandoDatabase_O,"CandoDatabase",core::CxxObject_O);
public:
  static CandoDatabase_sp make(core::Symbol_sp name);
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
  
  void	initialize();

public:
  typedef gc::SmallMap<core::Symbol_sp,FrameRecognizer_sp>	FrameRecognizers;
  typedef gc::SmallMap<core::Symbol_sp, Entity_sp>		Entities;
private: // archive
  core::Symbol_sp				_Name;
  string				_DateCreated;
  string				_DateUpdated;
  core::HashTable_sp          _Topologys;
  
private:
	/*! Return true if the entity with the name (name) is recognized and is a subclass of classId
	 */
  bool recognizesEntitySubClassOf(core::Symbol_sp name, core::Instance_sp mc);

	/*! Return true if the entity with the name (name) is recognized and of classId
	 */
  bool recognizesEntityOfClass(core::Symbol_sp name, core::Instance_sp mc);

	/*! Return a Cons of all entities that are of class (classId) */
  core::List_sp entitiesSubClassOfAsList(core::Instance_sp mc);


	/*! Return the entity with the given name and classId and if it doesn't exist throw an exception */
  Entity_sp getEntitySubClassOf(core::Symbol_sp name, core::Instance_sp mc);
 

	/*! Return the entity with the given name and classId */
  Entity_sp getEntityOfClass(core::Symbol_sp name, core::Instance_sp mc );


public:

  void	setDateCreated(const string& dc){this->_DateCreated=dc;};
  void	setDateUpdated(const string& dc){this->_DateUpdated=dc;};


  CL_LISPIFY_NAME("getName");
  CL_DEFMETHOD     core::Symbol_sp getName() { return this->_Name; };
  void	setName(core::Symbol_sp n) { this->_Name = n; };



	/*! Access contents using resource/name */
#ifdef USE_TOPOLOGY
  core::T_sp	oGetResource(core::Symbol_sp resource, core::Symbol_sp name );
#endif

  bool recognizesMonomerName(core::Symbol_sp nm);


	/*! Add a RepresentedEntityNameSet or Constitution */
  Entity_sp	addEntity(Entity_sp group);



  core::List_sp  representedEntityNameSetsAsList();

  core::List_sp	constitutionsAsList();


  bool recognizesNameOrPdb(core::Symbol_sp name);

//	Constitution_sp	constitutionForNameOrPdb(core::Symbol_sp name);
//	core::Symbol_sp constitutionNameForNameOrPdb(core::Symbol_sp name);
//	core::Symbol_sp getMonomerNameForNameOrPdb(core::Symbol_sp name);
//	core::Symbol_sp getPdbNameForNameOrPdb(core::Symbol_sp name);


  /*! Return (values topology t/nil) */
  core::T_mv findTopology(core::T_sp name,bool errorp) const;
  void setf_findTopology(core::T_sp name, core::T_sp topology);

  void walk_topologys(core::Function_sp func);
  

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------                                                                


public:


#ifdef USE_TOPOLOGY
  adapt::SymbolSet_sp getMonomersForSetOrConstitutionOrMonomerName(core::Symbol_sp nm);


  adapt::SymbolSet_sp expandEntityNameToTerminals(core::Symbol_sp nm) const;
  adapt::SymbolSet_sp expandEntityNamesToTerminals(adapt::SymbolSet_sp names) const;

  CL_LISPIFY_NAME("recognizesEntityName");
  CL_DEFMETHOD 	bool	recognizesEntityName(core::Symbol_sp nm) { return this->recognizesSetOrConstitutionOrMonomerName(nm);};


  Entity_sp getEntity(core::Symbol_sp nm) const;
#endif

#ifdef USE_TOPOLOGY
  adapt::SymbolSet_sp allUniqueCouplingNames();
#endif

		/*! When objects that depend on different CandoDatabase(s)
		 * are to be combined they have to reference a common CandoDatabase
		 * This function will signal to all of the objects that
		 * rely on this CandoDatabase that they need to move to another.
		 */
  void	giveYourDependantObjectsTo(CandoDatabase_sp newDb);


  core::List_sp	allMonomerNamesOrdered();
  void		testConsistency(std::ostream& out);
public:

  DEFAULT_CTOR_DTOR(CandoDatabase_O);
};


core::T_mv chem__findTopology(core::T_sp name, bool errorp = true);
}; // _Namespace chem

#endif
