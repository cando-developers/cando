     
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
#include "core/common.h"
#include "core/stringList.fwd.h"
#include "core/symbolSet.fwd.h"
#include "monomer.h"
#include "core/holder.h"
#include "addon/vector3.h"
#include "core/objRef.fwd.h"
#include "entity.fwd.h"
#include "chemPackage.h"


namespace chem {



SMART(Aggregate);
SMART(Entity);
SMART(SpecificContextSet);
SMART(Monomer);
SMART(ScorerVirtualMachine);
SMART(MultiMonomer);
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
class CandoDatabase_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,CandoDatabase_O,"CandoDatabase");
#if INIT_TO_FACTORIES
 public:
    static CandoDatabase_sp make(core::Symbol_sp name, core::Cons_sp frameRecognizers, core::Cons_sp representedEntityNameSets, core::Cons_sp constitutions);
#else
    DECLARE_INIT();
#endif
public:
//    void	archiveBase(core::ArchiveP node);
	void	initialize();

public:
	typedef	core::SymbolMap<FrameRecognizer_O>	FrameRecognizers;
	typedef	core::SymbolMap<MonomerCoordinates_O>		MonomerCoordinates;
	typedef core::SymbolMap<Entity_O>		Entities;
private: // archive
    core::Symbol_sp				_Name;
	string				_DateCreated;
	string				_DateUpdated;
	FrameRecognizers		_frameRecognizers;
	MonomerCoordinates		_MonomerCoordinates;
			/*! The Entities object maps
			 * names of Stereoisomers, Constitutions, RepresentedEntityNameSets
			 * and MonomerPacks to their respective objects.
			 * Any object that can be expanded ultimately into a list
			 * of Stereoisomer names is stored here.
			 */
	Entities	_Entities;
public:
	typedef MonomerCoordinates::iterator	monomerCoordinatesIterator;
private:
	/*! Return true if the entity with the name (name) is recognized and is a subclass of classId
	 */
	bool recognizesEntitySubClassOf(core::Symbol_sp name, core::Class_sp mc);

	/*! Return true if the entity with the name (name) is recognized and of classId
	 */
	bool recognizesEntityOfClass(core::Symbol_sp name, core::Class_sp mc);

	/*! Return a Cons of all entities that are of class (classId) */
	core::Cons_sp entitiesSubClassOfAsCons(core::Class_sp mc);


	/*! Return the entity with the given name and classId and if it doesn't exist throw an exception */
	Entity_sp getEntitySubClassOf(core::Symbol_sp name, core::Class_sp mc);
 

	/*! Return the entity with the given name and classId */
	Entity_sp getEntityOfClass(core::Symbol_sp name, core::Class_sp mc );


public:

	void	setDateCreated(const string& dc){this->_DateCreated=dc;};
	void	setDateUpdated(const string& dc){this->_DateUpdated=dc;};


    core::Symbol_sp getName() { return this->_Name; };
    void	setName(core::Symbol_sp n) { this->_Name = n; };



	/*! Access contents using resource/name */
    core::T_sp	oGetResource(core::Symbol_sp resource, core::Symbol_sp name );

	bool recognizesMonomerName(core::Symbol_sp nm);


	/*! Add a RepresentedEntityNameSet or Constitution */
	Entity_sp	addEntity(Entity_sp group);



	core::Cons_sp  representedEntityNameSetsAsCons();

	core::Cons_sp	constitutionsAsCons();


	bool recognizesNameOrPdb(core::Symbol_sp name);

	Constitution_sp	constitutionForNameOrPdb(core::Symbol_sp name);
	core::Symbol_sp constitutionNameForNameOrPdb(core::Symbol_sp name);
	core::Symbol_sp getMonomerNameForNameOrPdb(core::Symbol_sp name);
	core::Symbol_sp getPdbNameForNameOrPdb(core::Symbol_sp name);


	FrameRecognizer_sp	getFrameRecognizer(core::Symbol_sp nm)
    { _G();
	    FrameRecognizer_sp gr;
	    LOG(BF( "Agroup name in database=%d") % this->recognizesFrameRecognizerName(nm) );
            FrameRecognizers::iterator ig = this->_frameRecognizers.find(nm);
            if ( ig == this->_frameRecognizers.end() ) {
                SIMPLE_ERROR(BF("Could not find %s") % _rep_(nm));
            }
	    return ig->second;
	};


	void addFrameRecognizer(FrameRecognizer_sp rec);
	bool	recognizesFrameRecognizerName(core::Symbol_sp nm);


	bool	recognizesEntityNameSetName(core::Symbol_sp nm);

	bool	recognizesRepresentedEntityNameSet(core::Symbol_sp nm);



	
	/*! Return true if the Symbol references a MonomerPack */
	bool	recognizesMonomerPack(core::Symbol_sp nm);


	EntityNameSet_sp	getEntityNameSet(core::Symbol_sp nm);
	RepresentedEntityNameSet_sp getRepresentedEntityNameSet(core::Symbol_sp nm);


	bool	recognizesSetOrConstitutionOrMonomerName(core::Symbol_sp nm);


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------                                                                


public:


	core::SymbolSet_sp getMonomersForSetOrConstitutionOrMonomerName(core::Symbol_sp nm);


	core::SymbolSet_sp expandEntityNameToTerminals(core::Symbol_sp nm) const;
	core::SymbolSet_sp expandEntityNamesToTerminals(core::SymbolSet_sp names) const;

	bool	recognizesEntityName(core::Symbol_sp nm) { return this->recognizesSetOrConstitutionOrMonomerName(nm);};


	Entity_sp getEntity(core::Symbol_sp nm) const;


core::SymbolSet_sp allUniqueCouplingNames();

		/*! When objects that depend on different CandoDatabase(s)
		 * are to be combined they have to reference a common CandoDatabase
		 * This function will signal to all of the objects that
		 * rely on this CandoDatabase that they need to move to another.
		 */
	void	giveYourDependantObjectsTo(CandoDatabase_sp newDb);


    core::Cons_sp	allMonomerNamesOrdered();
    void		testConsistency(std::ostream& out);

public:

	monomerCoordinatesIterator begin_MonomerCoordinates_keyValue()
		{ return this->_MonomerCoordinates.begin(); };
	monomerCoordinatesIterator end_MonomerCoordinates_keyValue()
		{ return this->_MonomerCoordinates.end(); };

	monomerCoordinatesIterator	begin_MonomerCoordinates() {
			return this->_MonomerCoordinates.begin();
	};
	monomerCoordinatesIterator	end_MonomerCoordinates() {
			return this->_MonomerCoordinates.end();
	};

#if 0
	monomerCoordinatesValueIterator	begin_value_MonomerCoordinates() {
			return this->_MonomerCoordinates.begin_value();
	};
	monomerCoordinatesValueIterator	end_value_MonomerCoordinates() {
			return this->_MonomerCoordinates.end_value();
	};
#endif
		/*! Add the MonomerCoordinates but only
		 * with those context keys that this database requires
		 * Return the number of contexts that were added
		 */
	uint addMonomerCoordinates(MonomerCoordinates_sp coords);

	MonomerCoordinates_sp	getMonomerCoordinatesWithKey(MonomerName key);
	bool			recognizesMonomerCoordinatesKey(MonomerName key);

	core::Cons_sp uniqueMonomerCoordinatesAsCons();
	core::Cons_sp monomerCoordinatesKeysAsCons();

	SpecificContextSet_sp allSpecificMonomerContexts();

	MonomerCoordinates_sp	get(MonomerContext_sp context);
	bool recognizesContext(MonomerContext_sp context);

	void	saveAs(const string& fileName);

#if 0
	// I'll need to add this back once I figure out what Alchemists do
	void removeMonomerCoordinatesNotRequiredByAlchemists(core::Cons_sp alchemists);
#endif
public:

	DEFAULT_CTOR_DTOR(CandoDatabase_O);
};


}; // _Namespace chem

TRANSLATE(chem::CandoDatabase_O);
#endif
