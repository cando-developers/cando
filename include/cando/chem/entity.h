       
       
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
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,Entity_O,"Entity");
public:
	void initialize();
public:
//	void	archiveBase(core::ArchiveP node);
////	bool	loadFinalize(core::ArchiveP node);
private:
//	CandoDatabase_wp	_WeakCandoDatabase;
public:

	virtual string __repr__() const;

		/*! Return true if this name can not be expanded.
		 * For instance monomer groups can be expanded to every monomer within them.
		 */
    virtual bool	isTerminalName() { return false; };
    		/*! Get the name of this object
		 */
    virtual core::Symbol_sp	getName() const 	{_OF();SUBCLASS_MUST_IMPLEMENT();};

    		/*! Return a list of names that is expanded recursively as much as possible
		 */
    virtual adapt::SymbolSet_sp	expandedNameSet() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
    		/*! Return a list of names that is expanded recursively but terminate with
		 * representative names.
		 */
    /*! Recursively expand all of the EntityNames within us and accumulate a list of
      RepresentedEntityNameSets that have a representative and a set of Entities that
      that representative represents */
    virtual RepresentativeList_sp expandedRepresentativeList() const {_OF(); SUBCLASS_MUST_IMPLEMENT(); };

    /*! Call expandedRepresentativeList and then combine RepresentedEntityNameSets that
      have the same representative, this will give a minimal number of RepresentedEntityNameSets
    in terms of the number of representatives */
    RepresentativeList_sp minimalRepresentativeList() const;

		/*! If this object is a Stereoisomer then it has a Constitution otherwise it doesn't */
    virtual bool hasConstitution() { return false; };
    		/*! Return the constitution */
    virtual Constitution_sp constitution();


	Entity_O( const Entity_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(Entity_O);

};


};
TRANSLATE(chem::Entity_O);
#endif //]
