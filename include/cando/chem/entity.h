#ifdef USE_TOPOLOGY
/*
    File: entity.h
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
       
       
#ifndef	Entity_H //[
#define Entity_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/entity.fwd.h>
#include <cando/adapt/symbolSet.fwd.h>
//#include "geom/objectList.fwd.h"




#include <cando/chem/chemPackage.h>


namespace chem {



SMART(CandoDatabase);
SMART(Constitution);
SMART(Entity);
class Entity_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,Entity_O,"Entity",core::CxxObject_O);
 public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node) {}; // nothing to save
  
  void initialize();
 public:
 private:
//	CandoDatabase_sp	_WeakCandoDatabase;
 public:

  virtual string __repr__() const;

		/*! Return true if this name can not be expanded.
		 * For instance monomer groups can be expanded to every monomer within them.
		 */
  virtual bool	isTerminalName() { return false; };
    		/*! Return a list of names that is expanded recursively as much as possible
		 */
  virtual adapt::SymbolSet_sp	expandedNameSet() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
    		/*! Return a list of names that is expanded recursively but terminate with
		 * representative names.
		 */
    /*! Recursively expand all of the EntityNames within us and accumulate a list of
      RepresentedEntityNameSets that have a representative and a set of Entities that
      that representative represents */
  CL_LISPIFY_NAME("expandedRepresentativeList");
  CL_DEFMETHOD     virtual RepresentativeList_sp expandedRepresentativeList() const {_OF(); SUBCLASS_MUST_IMPLEMENT(); };

    /*! Call expandedRepresentativeList and then combine RepresentedEntityNameSets that
      have the same representative, this will give a minimal number of RepresentedEntityNameSets
    in terms of the number of representatives */
  RepresentativeList_sp minimalRepresentativeList() const;

		/*! If this object is a Stereoisomer then it has a Constitution otherwise it doesn't */
  CL_LISPIFY_NAME("hasConstitution");
  CL_DEFMETHOD     virtual bool hasConstitution() { return false; };
    		/*! Return the constitution */
  virtual Constitution_sp constitution();


  Entity_O( const Entity_O& ss ); //!< Copy constructor


  DEFAULT_CTOR_DTOR(Entity_O);

};


};
TRANSLATE(chem::Entity_O);
#endif //]
#endif
