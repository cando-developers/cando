/*
    File: monomerPack.cc
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
       
#define	DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/monomerPack.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/adapt/adapters.h>
#include <cando/adapt/symbolList.h>
#include <clasp/core/environment.h>
#include <cando/chem/candoDatabase.h>
#include <clasp/core/numbers.h>
#include <cando/adapt/stringList.h>
#include <cando/chem/alias.h>
#include <clasp/core/wrappers.h>

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


    
    
#define ARGS_chem__define_monomer_pack "(packName parts &optional atomAliases)"
#define DECL_chem__define_monomer_pack ""
#define DOCS_chem__define_monomer_pack "defineMonomerPack"
CL_LAMBDA(pack-name parts &optional atom-aliases);
CL_DEFUN core::T_sp chem__define_monomer_pack(core::Symbol_sp packName, core::List_sp parts, core::Cons_sp atomAliases )
    {
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


    
    
#define ARGS_chem__extend_aliases "(packName parts atomAliases)"
#define DECL_chem__extend_aliases ""
#define DOCS_chem__extend_aliases "extendAliases"
CL_DEFUN core::T_sp chem__extend_aliases(core::Symbol_sp packName, core::List_sp parts, core::Cons_sp atomAliases)
    {
	CandoDatabase_sp	bdb;
	core::T_sp		oarg1;
	bdb = getCandoDatabase();
	ASSERTNOTNULL(bdb);
	if ( !bdb->recognizesMonomerPack(packName) )
          SIMPLE_ERROR(BF("monomerPack("+_rep_(packName)+") has not been defined"));
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


    
    
#define ARGS_chem__set_monomer_pack "(packName parts &optional atomAliases)"
#define DECL_chem__set_monomer_pack ""
#define DOCS_chem__set_monomer_pack "setMonomerPack"
CL_LAMBDA(pack-name parts &optional atom-aliases);
CL_DEFUN core::T_sp chem__set_monomer_pack(core::Symbol_sp packName, core::List_sp parts, core::Cons_sp atomNames )
    {
	core::T_sp opack = chem__define_monomer_pack(packName,parts,atomNames);
	MonomerPack_sp pack = (opack).as<MonomerPack_O>();
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
    this->_InterestingAtomAliases = adapt::SymbolList_O::create();
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



void	MonomerPack_O::defineContentsFromCons(core::List_sp atomAliases, core::List_sp parts)
{
  core::List_sp		p;
  adapt::SymbolList_sp	aliases;
  core::Symbol_sp	name;
	// Create a string list of all the atom aliases
  aliases = adapt::SymbolList_O::create();
  for ( auto p : atomAliases ) {
    name = oCar(p).as<core::Symbol_O>();
    aliases->append(name);
  }
  this->setInterestingAtomAliasesFromSymbolList(aliases);
  for ( auto p : parts ) {
    core::List_sp entry = oCar(p);
    if ( core::cl__length(entry) <1 ) SIMPLE_ERROR(BF("monomerPack Entry contains no monomer name"));
    core::Symbol_sp monomerName = oCar(entry).as<core::Symbol_O>();
    this->addMonomerName(monomerName);
    if ( aliases->size() > 0 )
    {
      if ( core::cl__length(entry) < 2 ) SIMPLE_ERROR(BF("You defined atom aliases to each monomer must have interesting atom names"));
      core::List_sp interestingAtomNames = oSecond(entry);
      this->setInterestingAtomNamesForMonomerNameFromCons(monomerName,interestingAtomNames);
    }
  }
}




void MonomerPack_O::extendAliases( core::List_sp atomAliases, core::List_sp parts)
{_OF();
    adapt::SymbolSet_sp extendMonomers = adapt::SymbolSet_O::create();
    uint numberOfAtomAliases = core::cl__length(atomAliases);
    for ( auto aai : atomAliases ) {
      this->_InterestingAtomAliases = core::Cons_O::create(aai,this->_InterestingAtomAliases);
    }
    for ( auto cur : parts ) {
      core::List_sp oneExtend = oCar(cur);
      if ( core::cl__length(oneExtend) != 2 )
	{
          SIMPLE_ERROR(BF("Each extendAliases entry must have two elements: %s") % _rep_(oneExtend) );
	}
	core::Symbol_sp monomerName = oCar(oneExtend).as<core::Symbol_O>();
	core::List_sp aliasAtoms = oCadr(oneExtend);
	extendMonomers->insert(monomerName);
	if ( !this->_AtomIndexers->recognizesMonomerName(monomerName) )
	{
          SIMPLE_ERROR(BF("MonomerPack("+_rep_(this->getName())+") does not recognize monomer name("+_rep_(monomerName)+")"));
	}
	AtomIndexer_sp indexer = this->_AtomIndexers->getAtomIndexerForMonomerName(monomerName);
	if ( core::cl__length(aliasAtoms) != numberOfAtomAliases )
	{
	    stringstream ss;
	    ss << "Mismatch between number of alias atoms ";
	    ss << _rep_(aliasAtoms);
	    ss << " for monomer(" << _rep_(monomerName);
	    ss << ") and the number of atom aliases ";
	    ss << _rep_(atomAliases) ;
	    SIMPLE_ERROR(BF("%s")%ss.str());
	}
	indexer->appendConsOfTexts(aliasAtoms);
    }
    if ( !extendMonomers->equal(this->_EntityNames) )
    {
	stringstream ss;
	ss << "There is a mismatch between the monomers in this extendAliases command in the MonomerPack(" << this->getName() << ")" << std::endl;
	ss << "extendAliases monomer names: " << _rep_(extendMonomers) <<std::endl;
	ss << "MonomerPack names: " << _rep_(this->_EntityNames);
	SIMPLE_ERROR(BF("%s")%ss.str());
    }
}




void	MonomerPack_O::setInterestingAtomNamesForMonomerName(core::Symbol_sp monomerName, const string& atomIndexerNames)
{ 
CandoDatabase_sp		bdb;
AtomIndexer_sp			indexer;
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
    core::Symbol_sp monomerName, core::List_sp names )
{ 
CandoDatabase_sp		bdb;
AtomIndexer_sp			indexer;
    bdb = getCandoDatabase();
    ASSERTNOTNULL(this->_AtomIndexers);
    if ( !this->recognizesMonomerName(monomerName) )
    {
        SIMPLE_ERROR(BF("Unrecognized monomer name: %s %s") % monomerName % this->description() );
    }
    indexer = AtomIndexer_O::create();
    indexer->setFromList(names);
    this->_AtomIndexers->set(monomerName,indexer);
    this->_checkAtomIndexers();
}

void	MonomerPack_O::setInterestingAtomNamesForMonomerNameFromCons(
			core::Symbol_sp monomerName, core::List_sp names )
{ 
CandoDatabase_sp		bdb;
AtomIndexer_sp			indexer;
    bdb = getCandoDatabase();
    ASSERTNOTNULL(this->_AtomIndexers);
    if ( !this->recognizesMonomerName(monomerName) )
    {
        SIMPLE_ERROR(BF("Unrecognized monomer name: %s %s") % monomerName % this->sharedThis<MonomerPack_O>()->description() );
    }
    indexer = AtomIndexer_O::create();
    for ( auto p : names ) {
      core::Symbol_sp	interest = oCar(p).as<core::Symbol_O>();
      indexer->addAtomName(interest);
    }
    this->_AtomIndexers->set(monomerName,indexer);
    this->_checkAtomIndexers();
}

void	MonomerPack_O::setMonomerNameOrPdb(core::Symbol_sp nm)
{
    this->EntityNameSet_O::setMonomerNameOrPdb(nm);
    this->setInterestingAtomNamesForMonomerName(nm,"");
}


void	MonomerPack_O::addMonomerName(core::Symbol_sp nm)
{
    this->EntityNameSet_O::addMonomerName(nm);
    this->setInterestingAtomNamesForMonomerName(nm,"");
}

CL_LISPIFY_NAME("removeMonomerName");
CL_DEFMETHOD void	MonomerPack_O::removeMonomerName(core::Symbol_sp nm)
{
    this->_removeMonomerName(nm);
    this->_AtomIndexers->eraseEntryForMonomer(nm);
}

string	MonomerPack_O::getInterestingAtomNamesForMonomerName(core::Symbol_sp nm)
{
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


void	MonomerPack_O::setInterestingAtomAliasesFromSymbolList(adapt::SymbolList_sp names)
{
  this->_InterestingAtomAliases = _Nil<core::T_O>();
  for ( auto ni = names->begin(); ni!=names->end(); ++ni ) {
    this->_InterestingAtomAliases = core::Cons_O::create(*ni,this->_InterestingAtomAliases);
  }
    LOG(BF("For MonomerPack(%s) setting aliases(%s)") % this->getName().c_str() % this->_InterestingAtomAliases->asString().c_str()  );
}



string MonomerPack_O::getInterestingAtomAliasesAsString()
{
  return _rep_(this->_InterestingAtomAliases);
}


AtomIndexer_sp MonomerPack_O::getAtomIndexerForMonomerName(core::Symbol_sp monomerName)
{
    ASSERT(this->_AtomIndexers->recognizesMonomerName(monomerName));
    return this->_AtomIndexers->getAtomIndexerForMonomerName(monomerName);
}








void	MonomerPack_O::_checkAtomIndexers()
{
    return; // do nothing for now
MapOfMonomerNamesToAtomIndexers_O::smart_ptr_type	ais;
AtomIndexer_O::smart_ptr_type			ai;
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
{
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
{
  for ( auto ai : this->_InterestingAtomAliases ) {
    if ( core::cl__equal(ai,alias->getAtomAlias()) ) return true;
  }
  return false;
}


int	MonomerPack_O::getInterestingAtomAliasIndex(Alias_sp alias)
{
 int i = 0;
 for ( auto ai : this->_InterestingAtomAliases ) {
   if ( core::cl__equal(ai,alias->getAtomAlias()) ) return i;
   ++i;
 }
 return -1;
}








};
