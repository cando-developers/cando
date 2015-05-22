       
       
//
// (C) 2004 Christian E. Schafmeister
//


#ifndef MonomerRecognizer_H
#define MonomerRecognizer_H
#include "core/common.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "core/symbolSet.fwd.h"
#include "bond.h"
#include "addon/vector3.h"
//#include	"conformation.h"
#include "atom.h"
#include "core/stringSet.fwd.h"
#include "core/stringList.fwd.h"
#include "core/symbol.fwd.h"
#include "core/symbolList.fwd.h"
#include "entity.h"
#include "statusTracker.h"


#include "chemPackage.h"

#include "core/symbolList.fwd.h"// entityNameSet.h wants SymbolList needs symbolList.fwd.h

namespace chem {

SMART(Alias);
SMART(AtomIndexer);
SMART(CandoDatabase);
SMART(SymbolList);

SMART(EntityNameSetBase);
class EntityNameSetBase_O : public Entity_O
{
    LISP_BASE1(Entity_O);
    LISP_CLASS(chem,ChemPkg,EntityNameSetBase_O,"EntityNameSetBase");
#if INIT_TO_FACTORIES
 public:
    static EntityNameSetBase_sp make(core::Cons_sp entityNamesCons);
#else
public:
    DECLARE_INIT();
#endif
public:
    void initialize();
    string description() const;
    friend class CandoDatabase_O;
public:
//    void	archiveBase(core::ArchiveP node);
public:
    typedef	core::SymbolSet_sp	EntityNames;
protected:
    /*!Ring closing plugs have optional neighbors */
    bool		_Optional;
    core::SymbolSet_sp	_EntityNames;
protected:
    void _clear();
    void _removeMonomerName(core::Symbol_sp nm);
public:
    static EntityNameSetBase_sp	create2(core::Symbol_sp nm, core::Lisp_sp lisp);
public:
//	virtual void setCandoDatabase(CandoDatabase_sp db);

    virtual core::SymbolSet_sp	expandedNameSet();
    virtual RepresentativeList_sp expandedRepresentativeList() const;

    /*! Expand other and merge its names with mine */
    void mergeEntityNames(Entity_sp other);

public: // Lisp creation

public:
    virtual string getMonomerNameWithAtoms(core::Symbol_sp nm) { return "--- no atom names for Entity_O ---"; };

    void setOptional(bool b) { this->_Optional = b; };
    bool getOptional() { return this->_Optional; };

    void	setGroupNames(core::SymbolSet_sp s);


    /*! Return all other monomer names that are not recognized by
     * this EntityNameSetBase
     */
    core::SymbolSet_sp	getUnrecognizedMonomerNames();

    /*!Set this EntityNameSetBase up containing one monomer and
     * give it the same name as the monomer
     */
    virtual void	setMonomerNameOrPdb(core::Symbol_sp nm);


    /*! Add the group name and immediatly expand it to monomer names
     */
    virtual void	addMonomerName(core::Symbol_sp nm );
    virtual void	addEntityName(core::Symbol_sp nm ) { this->addMonomerName(nm);};
    virtual void	addMonomerNameNoChecks(core::Symbol_sp nm );


    bool	recognizesMonomerName(core::Symbol_sp name);
    //! Converts nameOrPdb to monomer name and checks if we recognize it
    bool	recognizesNameOrPdb(core::Symbol_sp nameOrPdb );

    core::SymbolSet_sp	getUnExpandedNames();

    core::SymbolSet_sp	getMonomerNames();
    string		getMonomerNamesAsString();

    //! String together group and monomer names
    virtual string		getKey();


    /*! Return a StringSet that combines every equivalentName
      from this EntityNameSetBase with every
      equivalentName from the passed EntityNameSetBase.
      This generates a list of monomer pairs that can be
      coupled by a certain coupling.
    */
//    core::StringSet_sp combineMonomerNames(string ins, EntityNameSetBase_sp e);



    core::Symbol_sp	getOnlyMonomerName();

    virtual bool isFinalized() { return true; };

    /*! Some subclasses support InterestingAtomAliases
     * This method can be called to figure out which ones do
     */
    virtual	bool	supportsInterestingAtomAliases() { return false;};
    virtual void setInterestingAtomNamesForMonomerName(core::Symbol_sp monomerName, const string& atomIndexerNames) {_OF(); SUBCLASS_MUST_IMPLEMENT();};
    virtual string getInterestingAtomNamesForMonomerName(core::Symbol_sp nm) { return "";};
    virtual bool hasInterestingAtomAlias(Alias_sp alias);
    virtual int getInterestingAtomAliasIndex(Alias_sp alias);
    virtual AtomIndexer_sp getAtomIndexerForMonomerName(core::Symbol_sp nm) {_OF();SUBCLASS_MUST_IMPLEMENT();};
    virtual core::SymbolList_sp getInterestingAtomAliases();

    virtual string getInterestingAtomAliasesAsString() { return "";};

    core::StringList_sp	getMonomerNamesOrdered();
    core::StringList_sp	getUnrecognizedMonomerNamesOrdered();



    /*! Expand my EntityNames to terminal EntityNames */
    void expandToTerminalEntityNames();


    /*! Contract the Entity names */
    void contractEntityNames(core::SymbolSet_sp entityNameSets);


    string	testEntityNameSetBase() { return "testEntityNameSetBase"; };

    EntityNameSetBase_O( const EntityNameSetBase_O& emr );

    DEFAULT_CTOR_DTOR(EntityNameSetBase_O);
};




/*! All EntityNameSet adds to EntityNameSetBase is that it stores a name for the EntityNameSet
 */
SMART(EntityNameSet);
class EntityNameSet_O : public EntityNameSetBase_O
{
    friend class CandoDatabase_O;
    LISP_BASE1(EntityNameSetBase_O);
    LISP_CLASS(chem,ChemPkg,EntityNameSet_O,"EntityNameSet");
#if INIT_TO_FACTORIES
 public:
    static EntityNameSet_sp make(core::Symbol_sp name);
#else
    DECLARE_INIT();
#endif
public:
	void initialize();
	string description() const;
public:
	static EntityNameSet_sp	create2(core::Symbol_sp nm,core::Lisp_sp lisp);
public:
//    void	archiveBase(core::ArchiveP node);
protected:
    core::Symbol_sp		_Name;
public:
    void setName(core::Symbol_sp nm) { this->_Name = nm;};
	core::Symbol_sp getName() const { return this->_Name;};
	EntityNameSet_O( const EntityNameSet_O& emr );

	DEFAULT_CTOR_DTOR(EntityNameSet_O);
};



SMART(EntityNameSetWithCap);
class EntityNameSetWithCap_O : public EntityNameSetBase_O
{
    LISP_BASE1(EntityNameSetBase_O);
    LISP_CLASS(chem,ChemPkg,EntityNameSetWithCap_O,"EntityNameSetWithCap");
public:
    void	initialize();
    string description() const;
public:
//    void	archiveBase(core::ArchiveP node);
public:
//    static EntityNameSetWithCap_sp create(core::Lisp_sp e,CandoDatabase_sp db);
public:
//	virtual void setCandoDatabase(CandoDatabase_sp db);
    /*! Add the group name and immediately expand it to
     * a list of monomers */
    void	addGroupName(core::Symbol_sp nm );
    bool	recognizesGroupName(core::Symbol_sp nm);


    EntityNameSetWithCap_O( const EntityNameSetWithCap_O& emr );
    DEFAULT_CTOR_DTOR(EntityNameSetWithCap_O);
};




};










TRANSLATE(chem::EntityNameSetBase_O);
TRANSLATE(chem::EntityNameSet_O);
TRANSLATE(chem::EntityNameSetWithCap_O);
#endif


