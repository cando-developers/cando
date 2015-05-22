       
#define	DEBUG_LEVEL_NONE

#include "core/common.h"
#include "core/str.h"
#include "core/stringSet.h"
#include "monomerPack.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "core/symbolSet.h"
#include "core/symbolList.h"
#include "core/environment.h"
#include "candoDatabase.h"
#include "core/numbers.h"
#include "core/stringList.h"
#include "alias.h"
#include "core/wrappers.h"

namespace chem {







/*
__BEGIN_DOC( candoScript.monomerPack.createMonomerPack, createMonomerPack)
\scriptCmd{createMonomerPack}{name:text monomersAndInterestingAtomNames:list atomAliases:list }\par
\scriptCmd{createMonomerPack}{name:text monomersAndInterestingAtomNames:list}\par

Create a MonomerPack and put it into the database with the name: \scriptArg{name}. A MonomerPack is a group of Stereoisomers each of which has zero or more atom names associated with it that will be built by CANDO during rapid searching through sequence and conformational space. 

\begin{verbatim}
( aaGLu = [createMonomerPack "aaGlu"
    [parts
        [addPart 'glu(S) [aliasAtoms 'OE ] ]
        [addPart 'glu(R) [aliasAtoms 'OE ] ]
    ]
    [atomAliases 'carbO ]
] )
\end{verbatim}

The command names "parts" and "atomAliases" are aliases for the "list" command.
__END_DOC
 */


    
    
#define ARGS_af_defineMonomerPack "(packName parts &optional atomAliases)"
#define DECL_af_defineMonomerPack ""
#define DOCS_af_defineMonomerPack "defineMonomerPack"
    core::T_sp af_defineMonomerPack(core::Symbol_sp packName, core::Cons_sp parts, core::Cons_sp atomAliases )
    {_G();
	CandoDatabase_sp	bdb;
	core::T_sp		oarg1;
	bdb = getCandoDatabase();
#if 0
	if ( bdb->recognizesUserMonomerPack(packName->get()) )
	    SIMPLE_ERROR(BF("monomerPack("+packName->get()+") is already defined"));
#endif
	MonomerPack_sp monomerPack = MonomerPack_O::create();
	monomerPack->setName(packName);
	monomerPack->defineContentsFromCons(atomAliases, parts);
    	// 
	// Add the MonomerPack to the CandoDatabase
	//
	bdb->addEntity(monomerPack);
	return monomerPack;
    }



/*
__BEGIN_DOC( candoScript.monomerPack.extendAliases, extendAliases)
\scriptCmd{extendAliases}{Text::monomerPackName Cons::monomersAndInterestingAtomNames Cons::atomAliases}\par

Lookup a MonomerPack in the CandoDatabase and extend the interesting atom list.



\begin{verbatim}
extendAliases "allBis"
    (:
        (: 'glu(S) (aliasAtoms 'OE ) )
        (: 'glu(R) (aliasAtoms 'OE ) )
    )
    (atomAliases 'carbO )
) )
\end{verbatim}
__END_DOC
 */


    
    
#define ARGS_af_extendAliases "(packName parts atomAliases)"
#define DECL_af_extendAliases ""
#define DOCS_af_extendAliases "extendAliases"
    core::T_sp af_extendAliases(core::Symbol_sp packName, core::Cons_sp parts, core::Cons_sp atomAliases)
    {_G();
	CandoDatabase_sp	bdb;
	core::T_sp		oarg1;
	bdb = getCandoDatabase();
	ASSERTNOTNULL(bdb);
	if ( !bdb->recognizesMonomerPack(packName) )
	    SIMPLE_ERROR(BF("monomerPack("+packName->__repr__()+") has not been defined"));
	MonomerPack_sp monomerPack = bdb->getEntity(packName).as<MonomerPack_O>();
	monomerPack->extendAliases(atomAliases, parts);
	return _Nil<core::T_O>();
    }



/*
__BEGIN_DOC( candoScript.monomerPack.setMonomerPack, setMonomerPack)
\scriptCmd{setMonomerPack}{name:text monomersAndInterestingAtomNames:list atomAliases:list }\par
\scriptCmd{setMonomerPack}{name:text monomersAndInterestingAtomNames:list}\par

Create a MonomerPack and put it into the database with the name: \scriptArg{name}, also create a local variable with the name \sa{name} containing this MonomerPack. A MonomerPack is a group of Stereoisomers each of which has zero or more atom names associated with it that will be built by CANDO during rapid searching through sequence and conformational space. 

\begin{verbatim}
[setMonomerPack "aaGlu"
    [parts
        [addPart 'glu(S) [aliasAtoms 'OE ] ]
        [addPart 'glu(R) [aliasAtoms 'OE ] ]
    ]
    [atomAliases 'carbO ]
]
\end{verbatim}

The commands "parts" and "atomAliases" are aliases for the "list" command.
__END_DOC
 */


    
    
#define ARGS_af_setMonomerPack "(packName parts &optional atomAliases)"
#define DECL_af_setMonomerPack ""
#define DOCS_af_setMonomerPack "setMonomerPack"
    core::T_sp af_setMonomerPack(core::Symbol_sp packName, core::Cons_sp parts, core::Cons_sp atomNames )
    {_G();
	core::T_sp opack = af_defineMonomerPack(packName,parts,atomNames);
	MonomerPack_sp pack = downcast<MonomerPack_O>(opack);
	core::Symbol_sp sym = pack->getName();
	sym->defparameter(opack);
	return opack;
    }










//
// Constructor
//

void	MonomerPack_O::initialize()
{
    this->Base::initialize();
    this->_InterestingAtomAliases = core::SymbolList_O::create();
    this->_AtomIndexers = MapOfMonomerNamesToAtomIndexers_O::create();
}

//
// Destructor
//

#ifdef XML_ARCHIVE
    void	MonomerPack_O::archive(core::ArchiveP node)
{
    this->EntityNameSet_O::archive(node);
    node->attribute("interestingAtomAliases",this->_InterestingAtomAliases);
    node->attribute("atomIndexers",this->_AtomIndexers);
}
#endif



void	MonomerPack_O::defineContentsFromCons(core::Cons_sp atomAliases, core::Cons_sp parts)
{_G();
core::Cons_sp		p;
core::SymbolList_sp	aliases;
core::Symbol_sp	name;
	// Create a string list of all the atom aliases
    aliases = core::SymbolList_O::create();
    for ( p=atomAliases; p.notnilp(); p=p->cdr() )
    {
	name = p->car<core::Symbol_O>();
	aliases->append(name);
    }
    this->setInterestingAtomAliasesFromSymbolList(aliases);
    for ( p=parts; p.notnilp(); p = p->cdr() )
    {
	core::Cons_sp entry = p->car<core::Cons_O>();
	if ( entry->length() <1 ) SIMPLE_ERROR(BF("monomerPack Entry contains no monomer name"));
	core::Symbol_sp monomerName = entry->car<core::Symbol_O>();
	this->addMonomerName(monomerName);
	if ( aliases->size() > 0 )
	{
	    if ( entry->length() < 2 ) SIMPLE_ERROR(BF("You defined atom aliases to each monomer must have interesting atom names"));
	    core::Cons_sp interestingAtomNames = entry->listref<core::Cons_O>(1);
	    this->setInterestingAtomNamesForMonomerNameFromCons(monomerName,interestingAtomNames);
	}
    }
}




void MonomerPack_O::extendAliases( core::Cons_sp atomAliases, core::Cons_sp parts)
{_OF();
    core::SymbolSet_sp extendMonomers = core::SymbolSet_O::create();
    uint numberOfAtomAliases = atomAliases->length();
    this->_InterestingAtomAliases->appendConsOfStrings(atomAliases);
    for ( core::Cons_sp cur = parts; cur.notnilp(); cur=cur->cdr() )
    {
	core::Cons_sp oneExtend = cur->car<core::Cons_O>();
	if ( oneExtend->length() != 2 )
	{
	    SIMPLE_ERROR(BF("Each extendAliases entry must have two elements: "+oneExtend->__repr__() ));
	}
	core::Symbol_sp monomerName = oneExtend->car<core::Symbol_O>();
	core::Cons_sp aliasAtoms = oneExtend->ocadr().as<core::Cons_O>();
	extendMonomers->insert(monomerName);
	if ( !this->_AtomIndexers->recognizesMonomerName(monomerName) )
	{
	    SIMPLE_ERROR(BF("MonomerPack("+this->getName()->__repr__()+") does not recognize monomer name("+monomerName->__repr__()+")"));
	}
	AtomIndexer_sp indexer = this->_AtomIndexers->getAtomIndexerForMonomerName(monomerName);
	if ( aliasAtoms->length() != numberOfAtomAliases )
	{
	    stringstream ss;
	    ss << "Mismatch between number of alias atoms ";
	    ss << aliasAtoms->__repr__();
	    ss << " for monomer(" << monomerName->__repr__();
	    ss << ") and the number of atom aliases ";
	    ss << atomAliases->__repr__() ;
	    SIMPLE_ERROR(BF("%s")%ss.str());
	}
	indexer->appendConsOfTexts(aliasAtoms);
    }
    if ( !extendMonomers->equal(this->_EntityNames) )
    {
	stringstream ss;
	ss << "There is a mismatch between the monomers in this extendAliases command in the MonomerPack(" << this->getName() << ")" << std::endl;
	ss << "extendAliases monomer names: " << extendMonomers->__repr__() <<std::endl;
	ss << "MonomerPack names: " << this->_EntityNames->__repr__();
	SIMPLE_ERROR(BF("%s")%ss.str());
    }
}




void	MonomerPack_O::setInterestingAtomNamesForMonomerName(core::Symbol_sp monomerName, const string& atomIndexerNames)
{ _G();
core::StringSet_O::smart_ptr		monomerNames;
CandoDatabase_sp		bdb;
AtomIndexer_sp			indexer;
core::StringSet_O::iterator		it;
    bdb = getCandoDatabase();
    ASSERTNOTNULL(this->_AtomIndexers);
    if ( !this->recognizesMonomerName(monomerName) )
    {
        SIMPLE_ERROR(BF("Unrecognized monomer name: %s %s") % monomerName % this->description() );
    }
    indexer = AtomIndexer_O::create();
    indexer->setFromAtomNames(atomIndexerNames);
    this->_AtomIndexers->set(monomerName,indexer);
    this->_checkAtomIndexers();
}


void	MonomerPack_O::setInterestingAtomNamesForMonomerNameStringList(
			core::Symbol_sp monomerName, core::StringList_sp names )
{ _G();
core::StringSet_O::smart_ptr		monomerNames;
CandoDatabase_sp		bdb;
AtomIndexer_sp			indexer;
core::StringSet_O::iterator		it;
    bdb = getCandoDatabase();
    ASSERTNOTNULL(this->_AtomIndexers);
    if ( !this->recognizesMonomerName(monomerName) )
    {
        SIMPLE_ERROR(BF("Unrecognized monomer name: %s %s") % monomerName % this->description() );
    }
    indexer = AtomIndexer_O::create();
    indexer->setFromStringList(names);
    this->_AtomIndexers->set(monomerName,indexer);
    this->_checkAtomIndexers();
}

void	MonomerPack_O::setInterestingAtomNamesForMonomerNameFromCons(
			core::Symbol_sp monomerName, core::Cons_sp names )
{ _G();
core::StringSet_O::smart_ptr		monomerNames;
CandoDatabase_sp		bdb;
AtomIndexer_sp			indexer;
core::StringSet_O::iterator		it;
    bdb = getCandoDatabase();
    ASSERTNOTNULL(this->_AtomIndexers);
    if ( !this->recognizesMonomerName(monomerName) )
    {
        SIMPLE_ERROR(BF("Unrecognized monomer name: %s %s") % monomerName % this->sharedThis<MonomerPack_O>()->description() );
    }
    indexer = AtomIndexer_O::create();
    for ( core::Cons_sp p=names; p.notnilp(); p=p->cdr() )
    {
	core::Str_sp	interest = p->car<core::Str_O>();
	indexer->addAtomName(interest->get());
    }
    this->_AtomIndexers->set(monomerName,indexer);
    this->_checkAtomIndexers();
}

void	MonomerPack_O::setMonomerNameOrPdb(core::Symbol_sp nm)
{_G();
    this->EntityNameSet_O::setMonomerNameOrPdb(nm);
    this->setInterestingAtomNamesForMonomerName(nm,"");
}


void	MonomerPack_O::addMonomerName(core::Symbol_sp nm)
{_G();
    this->EntityNameSet_O::addMonomerName(nm);
    this->setInterestingAtomNamesForMonomerName(nm,"");
}

void	MonomerPack_O::removeMonomerName(core::Symbol_sp nm)
{_G();
    this->_removeMonomerName(nm);
    this->_AtomIndexers->eraseEntryForMonomer(nm);
}

string	MonomerPack_O::getInterestingAtomNamesForMonomerName(core::Symbol_sp nm)
{_G();
string	s;
    LOG(BF("Looking atoms for MonomerName(%s)") % nm.c_str() );
    if ( this->_AtomIndexers->recognizesMonomerName(nm) )
    {
        s = this->_AtomIndexers->getAtomIndexerForMonomerName(nm)->asString();
	LOG(BF("  Found atomNames(%s)") % s.c_str()  );
	return s;
    }
    return "";
}


void	MonomerPack_O::setInterestingAtomAliasesFromSymbolList(core::SymbolList_sp names)
{_G();
    this->_InterestingAtomAliases->clear();
    this->_InterestingAtomAliases->appendSymbolList(names);
    LOG(BF("For MonomerPack(%s) setting aliases(%s)") % this->getName().c_str() % this->_InterestingAtomAliases->asString().c_str()  );
}



string MonomerPack_O::getInterestingAtomAliasesAsString()
{_G();
    return this->_InterestingAtomAliases->asString();
}


AtomIndexer_sp MonomerPack_O::getAtomIndexerForMonomerName(core::Symbol_sp monomerName)
{_G();
    ASSERT(this->_AtomIndexers->recognizesMonomerName(monomerName));
    return this->_AtomIndexers->getAtomIndexerForMonomerName(monomerName);
}








void	MonomerPack_O::_checkAtomIndexers()
{_G();
    return; // do nothing for now
MapOfMonomerNamesToAtomIndexers_O::smart_ptr	ais;
AtomIndexer_O::smart_ptr			ai;
MapOfMonomerNamesToAtomIndexers_O::iterator	ii;
int						mostNumberOfAtoms;
    if ( this->_AtomIndexers.nilp() ) return;
    ais = this->_AtomIndexers;
    mostNumberOfAtoms = 0;
    for ( ii = ais->begin(); ii != ais->end(); ii++ )
    {
	if ( ii->second->numberOfAtomNames() > mostNumberOfAtoms )
	{
	    mostNumberOfAtoms = ii->second->numberOfAtomNames();
	}
    }
    if ( mostNumberOfAtoms == 0 )
    {
	LOG(BF("Resetting _AtomIndexers") );
	this->_AtomIndexers = _Nil<MapOfMonomerNamesToAtomIndexers_O>();
    }
}


string	MonomerPack_O::getMonomerNameWithAtoms(core::Symbol_sp nm)
{_G();
stringstream	ss;
AtomIndexer_sp	atomIndexer;
    ASSERTNOTNULL(this->_AtomIndexers);
    ss.str("");
    ss << nm;
    if ( !this->recognizesMonomerName(nm) )
    {
        SIMPLE_ERROR(BF("Unrecognized monomer name: %s %s")% nm % this->sharedThis<MonomerPack_O>()->description() );
    }
    if ( this->_AtomIndexers->recognizesMonomerName(nm) )
    {
        atomIndexer = this->_AtomIndexers->getAtomIndexerForMonomerName(nm);
        ss << " [" << atomIndexer->asString() << "]";
    }
    return ss.str();
}



bool	MonomerPack_O::hasInterestingAtomAlias(Alias_sp alias)
{_G();
    LOG(BF("MonomerPack: %s") % this->getName().c_str() );
    LOG(BF("Checking to see if it recognizes atom alias: %s") % alias->getAtomAlias().c_str()  );
    LOG(BF("The interesting atom aliases I recognize are(%s)") % this->_InterestingAtomAliases->asString().c_str()  );
    return this->_InterestingAtomAliases->contains(alias->getAtomAlias());
}


int	MonomerPack_O::getInterestingAtomAliasIndex(Alias_sp alias)
{_G();
core::StringList_O::iterator	si;
    return this->_InterestingAtomAliases->indexOf(alias->getAtomAlias());
}




void MonomerPack_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<MonomerPack_O>()
	.def("removeMonomerName",&MonomerPack_O::removeMonomerName)
//	.def("setInterestingAtomAliasesFromString",&MonomerPack_O::setInterestingAtomAliasesFromString)
    ;
    Defun(defineMonomerPack);
    Defun(extendAliases);
    Defun(setMonomerPack);
}

void MonomerPack_O::exposePython(core::Lisp_sp lisp)
{_G();
#ifdef	USEBOOSTPYTHON //[
//    boost::python::def("create_MonomerPack",&create_MonomerPack,_lisp);
    PYTHON_CLASS(ChemPkg,MonomerPack,"","",_lisp)
	.def("removeMonomerName",&MonomerPack_O::removeMonomerName)
//	.def("setInterestingAtomAliasesFromString",&MonomerPack_O::setInterestingAtomAliasesFromString)
    ;

#endif
}


EXPOSE_CLASS(chem,MonomerPack_O);

};


