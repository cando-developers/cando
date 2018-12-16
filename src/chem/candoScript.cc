/*
    File: candoScript.cc
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

#include <stdlib.h>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/ql.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/symbolTable.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/fileSystem.h>
//#include "core/xmlLoadArchive.h"
//#include "core/xmlSaveArchive.h"
#include <cando/adapt/stringList.h>
#include <clasp/core/designators.h>
#include <clasp/core/bundle.h>
#include <clasp/core/evaluator.h>
#include <cando/geom/coordinateArray.h>
#include <cando/chem/candoScript.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/monomer.h>
#include <cando/chem/oligomer.h>
#include <cando/chem/mol2.h>
#include <cando/chem/loop.h>
//#include "monomerCoordinates.h"
//#include "groupPart.h"
//#include "fragmentCoordinates.h"
#include <cando/chem/calculatePosition.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitution.h>
#include <cando/chem/chemPackage.h>
#include <clasp/core/wrappers.h>

#ifdef	READLINE
extern "C" char *readline( const char* prompt);
extern "C" void add_history(char* line);
#endif

namespace chem {


#ifdef XML_ARCHIVE
    void	OligomerPart_Base_O::archiveBase(core::ArchiveP node)
{
   // donothing
}
#endif

#ifdef XML_ARCHIVE
    void	OligomerPart_Monomer_O::archive(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("_MonomerId",this->_MonomerId);
    node->attribute("_GroupName",this->_GroupName);
    node->attributeIfNotNil("_MonomerAliases",this->_MonomerAliases);
    node->attributeIfNotDefault<string>("_Comment",this->_Comment,"");
}
#endif




/*!
 * Replace the alchemist database of this Lisp environment.
 * If there already was one then tell it to give all of its
 * dependent objects to the new CandoDatabase
 */
void setCandoDatabase(CandoDatabase_sp bdb, core::Lisp_sp lisp)
{
    LOG(BF("Setting *DATABASE* to database@%p") % bdb.get() );
    _sym_candoDatabase->defparameter(bdb);
}

#ifdef XML_ARCHIVE
/*!
 * Load a CandoDatabase into this lisp environment.
 */
    void	loadCandoDatabase(const string& fileName, uint verbosity, core::Lisp_sp lisp)
    {
	CandoDatabase_sp	bdb;
	core::XmlLoadArchive_sp	xml;
	xml = core::XmlLoadArchive_O::create();
	xml->setVerbosity(verbosity);
	xml->open(fileName);
	if ( xml->contains("candoDatabase") )
	{
	    core::T_sp obj = xml->get("candoDatabase");
	    if (!obj.isA<CandoDatabase_O>() )
	    {
		SIMPLE_ERROR(BF("The 'candoDatabase' object is not of class CandoDatabase"));
	    }
	    bdb = obj.as<CandoDatabase_O>();
	} else if ( xml->contains("only" ) )
	{
	    core::T_sp obj = xml->get("only");
	    if (!obj.isA<CandoDatabase_O>() )
	    {
		SIMPLE_ERROR(BF("The 'only' object is not of class CandoDatabase"));
	    }
	    bdb = obj.as<CandoDatabase_O>();
	} else
	{
	    SIMPLE_ERROR(BF("Could not find database in file: "+fileName ));
	}
	setCandoDatabase(bdb,lisp);
    }
#endif




Monomer_sp	OligomerPart_Monomer_O::createMonomer(CandoDatabase_sp bdb)
{
  Monomer_sp mon = Monomer_O::create();
  mon->setId(this->_MonomerId);
//  mon->setGroupName(this->_GroupName);
//  mon->setAliasesFromCons(this->_MonomerAliases);
  return mon;
}



#ifdef XML_ARCHIVE
void	OligomerPart_Link_O::archive(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("_Monomer1Id",this->_Monomer1Id);
    node->attribute("_Coupling",this->_Coupling);
    node->attribute("_Monomer2",this->_Monomer2);
}
#endif



Monomer_sp	OligomerPart_Link_O::createMonomer(CandoDatabase_sp bdb)
{
    return this->_Monomer2->createMonomer(bdb);
}



/*
__BEGIN_DOC( candoScript.database, section, Database commands)
These commands load and inspect the {\CANDOSCRIPT} environment database. There is only
one database loaded at a time.
__END_DOC
*/

/*
__BEGIN_DOC(candoScript.database.database,subsection,database)
\scriptCmd{database}{fileName:text}

Load the database from \scriptArg{fileName}.
__END_DOC
*/
#define ARGS_chem__database "(fileName &optional (verbosity 0))"
#define DECL_chem__database ""
#define DOCS_chem__database "database"
CL_LAMBDA(fileName &optional (verbosity 0));
CL_DEFUN core::T_sp chem__database(core::T_sp fileNameDesig, core::Fixnum_sp overbosity)
    {
	IMPLEMENT_ME();
#ifdef XML_ARCHIVE
	core::Path_sp path = core::coerce::pathDesignator(fileNameDesig);
	uint verbosity = overbosity->get();
	LOG(BF("Opening file(%s)") % path->asString() );
	loadCandoDatabase(path->asString(),verbosity);
#endif
	return _Nil<core::T_O>();
    }

/*
__BEGIN_DOC(candoScript.database.standardDatabase,standardDatabase)
\args{((Text databaseFileName) &optional (Bool loadSeed) (Int (verbosity 0)))}
\returns{()}

Load the database with the name \emph{databaseFileName}.
If \emph{loadSeed} is true then it loads the seed database (no MonomerCoordinates) 
- otherwise it loads the standard database.
If :verbosity is set to a non zero value then debugging information is printed during loading.
__END_DOC
*/

    
    
#define ARGS_chem__standard_database "(pathDesig &optional loadSeed (overbosity 0))"
#define DECL_chem__standard_database ""
#define DOCS_chem__standard_database "standardDatabase"
CL_LAMBDA(pathDesig &optional loadSeed (overbosity 0));
CL_DEFUN core::T_sp chem__standard_database(core::T_sp pathDesig, core::T_sp loadSeed, core::Fixnum_sp overbosity)
    {
	IMPLEMENT_ME();
#ifdef XML_ARCHIVE
	core::Path_sp path = core::coerce::pathDesignator(pathDesig);
	uint verbosity = overbosity->get();
	stringstream name;
	if ( loadSeed.isTrue() )
	{
	    name << "seed_" << path->asString() << ".cxml";
	} else
	{
	    name << "database_" << path->asString() << ".cxml";
	}
	boost_filesystem::path filePath  = _lisp->bundle().getDatabasesDir() / name.str();
	LOG(BF("About to load alchemist database(%s)") % filePath.string());
	loadCandoDatabase(filePath.string(),verbosity);
#endif
	return _Nil<core::T_O>();
    }

/*
__BEGIN_DOC(candoScript.database.bundleDatabasePath,bundleDatabasePath)
\scriptCmd{bundleDatabasePath}{directoryName:text}

Return the full path of a file in the bundle database directory.
__END_DOC
*/
#define ARGS_chem__bundle_database_path "(pathDesig)"
#define DECL_chem__bundle_database_path ""
#define DOCS_chem__bundle_database_path "bundleDatabasePath"
CL_DEFUN core::T_sp chem__bundle_database_path(core::T_sp pathDesig)
{
    IMPLEMENT_ME();
#ifdef XML_ARCHIVE
    core::Path_sp pathName = core::coerce::pathDesignator(pathDesig);
    boost_filesystem::path filePath  = _lisp->bundle().getDatabasesDir() / pathName->asString();
    core::Str_sp fullPathName = core::Str_O::create(filePath.string());
    return fullPathName;
#endif
    return _Nil<core::T_O>();
}


/*
__BEGIN_DOC(candoScript.database.setDatabase,setDatabase)
\scriptCmd{standardDatabase}{directoryName:text}

Set the database.
__END_DOC
*/
core::T_sp	prim_setDatabase(core::Function_sp e, core::Cons_sp args, core::Environment_sp environ, core::Lisp_sp lisp)
{
CandoDatabase_sp	cdb;
    cdb = args->car<CandoDatabase_O>();
    setCandoDatabase(cdb,lisp);
    return _Nil<core::T_O>();
}



#if 0
#if !KINEMATICS
/*
__BEGIN_DOC(candoScript.database.contextGrep)
\positional{Text::contextKeySubstring}

Search for contexts with keys that contain the substring.
__END_DOC
*/
core::T_sp	prim_contextGrep(core::Function_sp e, core::Cons_sp args, core::Environment_sp environ, core::Lisp_sp lisp)
{
string	part;
    CandoDatabase_sp bdb = getCandoDatabase();
    if ( bdb.nilp() )
    {
	SIMPLE_ERROR(BF("You must define a AlchemistDatabse"));
    }
    part = args->car<core::Str_O>()->get();

    CandoDatabase_O::monomerCoordinatesIterator mit;
    for ( mit=bdb->begin_MonomerCoordinates_keyValue();
    		mit != bdb->end_MonomerCoordinates_keyValue(); mit++ )
    {
	MonomerContext_sp context = mit->second->getContext();
	adapt::StringSet_sp allKeys = context->getAllSpecificKeys();
	uint times = 0;
	string firstFind;
        adapt::StringSet_O::iterator sit;
	for ( sit=allKeys->begin(); sit!=allKeys->end(); sit++ )
	{
	    if ( (*sit).find(part) != string::npos )
	    {
		firstFind = (*sit);
		times++;
	    }
	}
	if ( times > 0 )
	{
	    lisp->print(BF( "%s") % mit->second->getComment().c_str() );
	    lisp->print(BF("     Found %d times first: %s") %times % firstFind.c_str() );
	} 
    }
    return _Nil<core::T_O>();
}
#endif
#endif



#if 0
/*
__BEGIN_DOC(candoScript.RequiredPlug.RequiredPlug,subsection,RequiredPlug)
\scriptCmd{RequiredPlug}{plug: plug}

Return the plug. This is a function that passes through the Plug when we aren't useing RequiredPlug objects but rather directly using Plug objects.
__END_DOC
*/
#define ARGS_prim_RequiredPlug "(plug)"
#define DECL_prim_RequiredPlug ""
#define DOCS_prim_RequiredPlug ""
core::T_sp	prim_RequiredPlug(core::Function_sp e, core::Cons_sp args, core::Environment_sp environ, core::Lisp_sp lisp)
{
    Plug_sp p = environ->lookup(ChemPkg,"plug").as<Plug_O>();
    return p;
}

/*
__BEGIN_DOC(candoScript.database.describeDatabase,subsection,describeDatabase)
\scriptCmd{describeDatabase}{}

Describe the contents of the database.
__END_DOC
*/
core::T_sp	prim_describeDatabase(core::Function_sp e, core::List_sp args, core::Environment_sp environ, core::Lisp_sp lisp)
{
core::Str_sp		fileName;
CandoDatabase_sp	bdb;
    bdb = getCandoDatabase();
    core::List_sp cur = bdb->constitutionsAsList();
    for ( ; cur.notnilp(); cur = cur->cdr() )
    {
	Constitution_sp con = core::oCar(cur).as<Constitution_O>();
	string stereoisomers = con->getMonomerNamesAsStringList()->asString();
	lisp->print(BF(" Constitution: %12s -- %s") % _rep_(con->getName()) % stereoisomers );
    }
    return _Nil<core::T_O>();
}
#endif



/*
__BEGIN_DOC(candoScript.database.createVirtualAtom,subsection,createVirtualAtom)
\scriptCmd{createVirtualAtom}{name:String::virtualAtomName) stereoIsomer:String::stereoIsomerName fragmentName:String::fragmentName position:'Cons::positionCalculator}

Create a virtual atom and insert it into the fragment with name \sa{fragmentName} for very MonomerCoordinate with stereoisomer name \sa{stereoIsomerName}.
If \sa{fragmentName} isn't specified then CANDO puts the virtual atom in the same fragment as the closest atom as defined by \sa{positionCalculator}.

A VirtualAtom is a named point in space calculated from other atoms in a MonomerCoordinates. They let us do things like precalculate the point where a metal atom would coordinate off of a nitrogen ligand in a Histidine relative to the imidazole ring. 
__END_DOC
*/



#define ARGS_chem__create_virtual_atom "()"
#define DECL_chem__create_virtual_atom ""
#define DOCS_chem__create_virtual_atom "createVirtualAtom"
CL_DEFUN core::T_sp chem__create_virtual_atom()
{
    IMPLEMENT_ME(); // Handle new MonomerCoordinates
#if !KINEMATICS
    core::KeyedArguments_sp kargs = core::KeyedArguments_O::createFromKeyedObjectCons(args,lisp);
    string virtualAtomName = kargs->getStringAndRemove("name");
    string constitutionName = kargs->getStringAndRemove("constitutionName");
    string fragmentName = kargs->getStringAndRemoveOrDefault("fragmentName","");
    CalculatePosition_sp calcPos= kargs->getAndRemove("position").as<CalculatePosition_O>();
    CandoDatabase_sp bdb = getCandoDatabase();
    CandoDatabase_O::monomerCoordinatesIterator mi;
    uint added = 0;
    uint collisions = 0;
    //
    // Add the virtual atom to every constitution with this monomer name
    //
    {_BLOCK_TRACEF(BF("Adding virtual atom(%s) to constitution(%s)") % virtualAtomName.c_str()%constitutionName.c_str() );
	CandoDatabase_O::constitutionIterator ci;
	for ( ci=bdb->begin_Constitutions();
	      ci!=bdb->end_Constitutions(); ci++ )
	{
	    LOG(BF("Looking at constitution(%s)") % (*ci)->getName() );
	    if ( (*ci)->getName() == constitutionName )
	    {
		LOG(BF("Found constitution with name(%s)") % constitutionName );
		ConstitutionAtoms_sp res = (*ci)->getConstitutionAtoms();
		ConstitutionVirtualAtom_sp virt = ConstitutionVirtualAtom_O::create(virtualAtomName,calcPos);
		res->addConstitutionVirtualAtom(virt);
		added++;
	    } else
	    {
		LOG(BF("Does not contain stereoisomer") );
	    }
	}
    }
    if ( added == 0 )
    {
	SIMPLE_ERROR(BF("Could not find constitution: "+constitutionName));
    }
    //
    // Add the virtual atom to all MonomerCoordinates
    //
    {_BLOCK_TRACEF(BF("Adding VirtualAtom(%s) to MonomerCoordinates") % virtualAtomName.c_str() );
        gctools::SmallOrderedSet<MonomerCoordinates_sp>	uniqueMonomerCoordinates;
	for ( mi=bdb->begin_MonomerCoordinates(); 
	      mi!=bdb->end_MonomerCoordinates();
	      mi++ )
	{ 
	    uniqueMonomerCoordinates.insert(mi->second);
	}
        gctools::SmallOrderedSet<MonomerCoordinates_sp>::iterator ui;
        gctools::SmallOrderedSet<MonomerCoordinates_sp>	interestingMonomerCoordinates;
	for ( ui=uniqueMonomerCoordinates.begin();
	      ui!=uniqueMonomerCoordinates.end();
	      ui++ )
	{ _BLOCK_TRACEF(BF("Checking MonomerCoordinates with context: %s") % (*ui)->getContext()->getKey().c_str() );
	    MonomerContext_sp context = (*ui)->getContext();
	    EntityNameSetBase_sp focus = context->getFocus();
	    string focusMonomerName = focus->getOnlyMonomerName();
	    CandoDatabase_sp bdb = getCandoDatabase();
	    Constitution_sp constitution = bdb->constitutionForNameOrPdb(focusMonomerName);
	    if ( constitution->getName() == constitutionName )
	    {
		LOG(BF( "The context focus monomer(%s) has the constitution(%s)")%
		    focus->getMonomerNamesAsString() % constitutionName );
		interestingMonomerCoordinates.insert(*ui);
	    }
	}
	for ( gctools::SmallOrderedSet<MonomerCoordinates_sp>::iterator si=interestingMonomerCoordinates.begin();
	      si!=interestingMonomerCoordinates.end(); si++ )
	{ _BLOCK_TRACEF(BF("Checking MonomerCoordinates with context: %s") % _rep_((*si)->getContext()) );
	    // If the fragmentName wasn't defined then use the
	    // CalculatePosition object to tell us what atom is closest
	    // to the virtual atom and put it in the same fragment
	    // as that
	    MonomerCoordinates_sp monCoords = *si;
	    FragmentCoordinates_sp fragmentCoords;
	    if ( fragmentName == "" )
	    {
		fragmentCoords = monCoords->getFragmentCoordinatesThatContainAtomWithName(calcPos->closestAtomName());
		LOG(BF( "Adding virtual atom(%s) to fragmentCoords for fragment(%s)")%
		    virtualAtomName % fragmentCoords->getFragmentName() );
		LOG(BF( "For the virtual atom(%s) the closestAtomName was(%s)")%
		    virtualAtomName % calcPos->closestAtomName() );
	    } else
	    {
		fragmentCoords = monCoords->getFragmentCoordinatesWithName(fragmentName);
		if ( fragmentCoords.nilp() )
		{
		    SIMPLE_ERROR(BF("In createVirtualAtom: Unrecognized fragmentCoordinates name("+fragmentName+")"));
		}
	    }
	    // Here create the Virtual atom
	    //
	    LOG(BF( "Creating virtual atom in FragmentCoordinates(%s)")% 
		fragmentCoords->getFragmentName() );
	    if ( !fragmentCoords->addVirtualAtom(virtualAtomName,calcPos,lisp) )
		collisions++;
	}
    }
    lisp->print(BF("Added virtual atom(%s) for constitution(%s)") % virtualAtomName % constitutionName );
    if ( collisions > 0 )
    {
        lisp->print(BF(" saw %d collisions") % collisions );
    }
    lisp->print(BF(""));
    return _Nil<core::T_O>();
#endif
}

/*
__BEGIN_DOC( candoScript.general) % section, General commands)
These commands apply to a variety of objects.
__END_DOC
*/

/*
__BEGIN_DOC( candoScript.general.saveWithCandoDatabase, subsection, saveWithCandoDatabase)
\scriptCmd{save}{Object::object Text::filename}

Save the \sa{object} along with the system CandoDatabase to the \sa{filename} in Cando-XML format.
__END_DOC
 */




#define ARGS_chem__save_archive_with_auto_set_cando_database "()"
#define DECL_chem__save_archive_with_auto_set_cando_database ""
#define DOCS_chem__save_archive_with_auto_set_cando_database "saveArchiveWithAutoSetCandoDatabase"
CL_DEFUN core::T_sp chem__save_archive_with_auto_set_cando_database()
{
    IMPLEMENT_ME();
#if 0 //fix below
    core::KeyedArguments_sp kargs = core::KeyedArguments_O::createFromKeyedObjectCons(args,lisp);
    if ( kargs->numberOfPositionalArguments()!=2 ) SIMPLE_ERROR(BF("You must provide the object and file name"));
    core::T_sp obj = kargs->getPositionalArgument(0).as<core::T_O>();
    core::Str_sp fileName = kargs->getPositionalArgument(1).as<core::Str_O>();
    int debugLevel = kargs->getIntAndRemoveOrDefault("verbosity",0);
    core::XmlSaveArchive_sp archive = core::XmlSaveArchive_O::create();
    archive->setVerbosity(debugLevel);
#if 0 // Archive can't do this anymore
    archive->putCandoDatabase(getCandoDatabase());
#endif
    archive->put("only",obj);
    archive->saveAs(fileName->get());
//    obj->saveXmlAs(fileName->get(),"only");
    return _Nil<core::T_O>();
#endif
}


/*
__BEGIN_DOC( candoScript.general.loadMol2, subsection, loadMol2 )
\scriptCmdRet{loadMol2}{Text::fileName}{Matter::result}

Loads the Tripos ``mol2'' file and returns the Aggregate within it.
__END_DOC
 */

CL_DEFUN geom::SimpleVectorCoordinate_sp chem__make_simple_vector_coordinate_from_atom_list(core::List_sp atoms)
{
  size_t num = core::cl__length(atoms);
  geom::SimpleVectorCoordinate_sp coords = geom::SimpleVectorCoordinate_O::make(num);
  int idx(0);
  for ( auto cur : atoms ) {
    Atom_sp a = gc::As<chem::Atom_sp>(oCar(cur));
    Vector3 v = a->getPosition();
    (*coords)[idx] = v;
    ++idx;
  }
  return coords;
}

CL_DOCSTRING(R"doc(Load the only or the first aggregate from the mol2 file.)doc");
CL_DEFUN core::T_sp chem__load_mol2(core::T_sp fileName)
{
  Mol2File fin;
  fin.openFileName(fileName);
  Aggregate_sp agg = mol2Read(fin);
  return agg;
}


CL_DOCSTRING(R"doc(Load all or a number of aggregates from the mol2 file.)doc");
CL_LAMBDA(file-name &optional number-to-load);
CL_DEFUN core::T_sp chem__load_mol2_list(core::T_sp fileName, core::T_sp number_to_load)
{
  core::T_sp make_progress;
  core::T_sp progress_advance;
  core::T_sp progress_done;
  if (chem_verbose(0)) {
    make_progress = core::eval::funcall(cl::_sym_findSymbol,core::SimpleBaseString_O::make("MAKE-PROGRESS-BAR"),
                                        core::SimpleBaseString_O::make("CANDO"));
    progress_advance = core::eval::funcall(cl::_sym_findSymbol,core::SimpleBaseString_O::make("PROGRESS-ADVANCE"),
                                           core::SimpleBaseString_O::make("CANDO"));
    progress_done = core::eval::funcall(cl::_sym_findSymbol,core::SimpleBaseString_O::make("PROGRESS-DONE"),
                                        core::SimpleBaseString_O::make("CANDO"));
    if (make_progress.nilp()||progress_advance.nilp()||progress_done.nilp()) {
      SIMPLE_ERROR(BF("Could not get progress bar functions make-progress, progress-advance, progress-done"));
    }
  }
  Mol2File fin;
  fin.openFileName(fileName);
  core::T_sp progress_bar = _Nil<core::T_O>();
  if (chem_verbose(0)) {
    if (number_to_load.nilp()) {
      progress_bar = core::eval::funcall(make_progress,
                                         INTERN_(kw,total), core::cl__file_length(fin.fIn));
    } else {
      progress_bar = core::eval::funcall(make_progress,
                                         INTERN_(kw,total), core::make_fixnum((size_t)number_to_load));
    }
  }
  ql::list result;
  size_t count = ~0;
  if (number_to_load.fixnump()) {
    count = number_to_load.unsafe_fixnum();
  }
  size_t num = 0;
  while (!fin.eof()) {
    if (count) {
      core::T_sp tagg = mol2Read(fin);
      if (gc::IsA<Aggregate_sp>(tagg)) {
        result << tagg;
        if (chem_verbose(0)) {
          if (number_to_load.nilp()) {
            core::eval::funcall(progress_advance, progress_bar, core::cl__file_position(fin.fIn,_Nil<core::T_O>()));
          } else {
            core::eval::funcall(progress_advance, progress_bar, core::make_fixnum((size_t)num));            
          }
        }
      }
    } else {
      break;
    }
    ++num;
    --count;
  }
  if (chem_verbose(0)) {
    core::eval::funcall(progress_done,progress_bar);
  }
  return result.cons();
}


/*
__BEGIN_DOC( candoScript.general.saveMol2, subsection, saveMol2)
\scriptCmd{saveMol2}{ Matter::structure Text::fileName}

Saves \sa{structure} to the Tripos "mol2" file with the name: \sa{fileName}.
__END_DOC
 */


CL_LAMBDA(matter dest-desig &optional use-sybyl-types);
CL_DEFUN core::T_sp chem__save_mol2(Matter_sp matter, core::T_sp destDesig, bool useSybylTypes)
{
  mol2WriteMatterToFileName(matter,destDesig,useSybylTypes);
  return _Nil<core::T_O>();
}


CL_LAMBDA(matter &optional use-sybyl-types);
CL_DEFUN std::string chem__aggregate_as_mol2_string(Aggregate_sp matter, bool useSybylTypes)
{
  stringstream ss;
  mol2WriteAggregateStream(matter,ss,useSybylTypes);
  return ss.str();
}

/*!
 * Look for a residue using an identifier,
 * this looks for residues with a sequence number if the identifier is an int
 * or a residue name or a residue alias
 */

Residue_sp findResidue(Matter_sp matter, core::T_sp residueIdentifier )
{
Molecule_sp molecule;
core::Fixnum_sp	residueSequenceNumber;

	    // downcast the identifer to a Symbol and an Int object
	    // one of them will be nil and the other will have a value
    core::Symbol_sp resIdName = residueIdentifier.asOrNull<core::Symbol_O>();
    core::Fixnum_sp resIdSeqNum = residueIdentifier.asOrNull<core::Fixnum_O>();
    Loop l;
    l.loopTopGoal(matter,RESIDUES);
    while ( l.advanceLoopAndProcess() )
    {
        Residue_sp res = l.getResidue();
        LOG(BF("Looking at residue with sequence number: %d") % res->getFileSequenceNumber()  );
	bool foundResidue = false;
	if ( resIdName )
	{
          LOG(BF("Checking if residue has name(%s) that matches(%s)") % res->getName().c_str() % _rep_(resIdName) );
	    if ( resIdName==res->getName() )
	    {
              LOG(BF("Found residue with sequence name: %s") % _rep_(resIdName) );
		foundResidue = true;
	    } else if ( res->recognizesMonomerAlias(resIdName) )
	    {
              LOG(BF("Found residue with Monomer alias: %s") % _rep_(resIdName) );
		foundResidue = true;
	    }
	} else if ( resIdSeqNum )
	{
	    LOG(BF("Checking if residue has fileSequenceNumber(%d) that matches(%d)") % res->getFileSequenceNumber() % resIdSeqNum->get()  );
	    if ( (int)(res->getFileSequenceNumber()) == resIdSeqNum.unsafe_fixnum())
	    {
	        LOG(BF("Found residue with sequence number: %s") % resIdSeqNum->get()  );
		foundResidue = true;
	    }
	}
	if ( foundResidue )
	{
	    return res;
	}
    }
    return _Nil<Residue_O>();
}

/*
__BEGIN_DOC( candoScript.general.findResidue, subsection, findResidue)
\scriptCmdRet{findResidue}{ structure:Molecule residueIdentifier}{Residue}\par
\scriptCmdRet{findResidue}{ structure:Aggregate chain:text residueIdentifier}{Residue}

Return the Residue in the Aggregate or Molecule: \scriptArg{structure}. The residueIdentifier can either be a sequence number of a residue name.
__END_DOC
 */




#define ARGS_chem__find_residue "(&rest args)"
#define DECL_chem__find_residue ""
#define DOCS_chem__find_residue "findResidue"
CL_LAMBDA(&rest args);
CL_DEFUN core::T_sp chem__find_residue(core::List_sp args)
{
    Molecule_sp molecule;
    core::Fixnum_sp	residueSequenceNumber;
    core::T_sp residueIdentifier;
    if ( core::cl__length(args)==3 ) 
    {
      Aggregate_sp agg = args.asCons()->onth(0).as<Aggregate_O>();
      core::String_sp chain = args.asCons()->onth(1).as<core::String_O>();
      residueIdentifier = args.asCons()->onth(2).as<core::T_O>();
      molecule = gc::As<Molecule_sp>(agg->contentWithName(chemkw_intern(chain->get_std_string())));
    } else if ( core::cl__length(args)==2 ) 
    {
      molecule = args.asCons()->onth(0).as<Molecule_O>();
      residueIdentifier = args.asCons()->onth(1).as<core::T_O>();
    } else
    {
    	SIMPLE_ERROR(BF("You must provide a molecule, residueId" ));
    }
    Residue_sp res = findResidue(molecule,residueIdentifier);
    if ( res.nilp() )
    {
	stringstream serr;
	serr << "Molecule does not contain residue with identifier: " << _rep_(residueIdentifier).c_str() ;
	SIMPLE_ERROR(BF(serr.str()));
    }
    return res;
}






/*
__BEGIN_DOC( candoScript.general.atomPos, subsection, atomPos)
\scriptCmdRet{atomPos}{ structure:Molecule residueNumber::int atomName:text }{vector}\par
\scriptCmdRet{atomPos}{ structure:Aggregate chain:text residueNumber::int atomName:text}{vector}

Return the position of the atom in the Aggregate or Molecule: \scriptArg{structure}. Identify the atom with \scriptArg{residueNumber} and \sa{atomName} and optionally the \sa{molecule/chain}. This function returns the position as a three element list of numbers.
__END_DOC
 */



#define ARGS_chem__atom_pos "(&rest args)"
#define DECL_chem__atom_pos ""
#define DOCS_chem__atom_pos "atomPos"
CL_LAMBDA(&rest args);
CL_DEFUN core::T_sp chem__atom_pos(core::List_sp args)
{
    Molecule_sp molecule;
    core::Fixnum_sp	residueSequenceNumber;
    core::Symbol_sp	atomName;
    core::T_sp residueIdentifier;
    if ( core::cl__length(args)==4 ) 
    {
      Aggregate_sp agg = args.asCons()->onth(0).as<Aggregate_O>();
      core::String_sp chain = args.asCons()->onth(1).as<core::String_O>();
      residueIdentifier = args.asCons()->onth(2).as<core::T_O>();
      atomName = args.asCons()->onth(3).as<core::Symbol_O>();
      molecule = gc::As<Molecule_sp>(agg->contentWithName(chemkw_intern(chain->get_std_string())));
    } else if ( core::cl__length(args)==3 ) 
    {
      molecule = args.asCons()->onth(0).as<Molecule_O>();
      residueIdentifier = args.asCons()->onth(1).as<core::T_O>();
      atomName = args.asCons()->onth(2).as<core::Symbol_O>();
    } else
    {
    	SIMPLE_ERROR(BF("You must provide a molecule, residueId and atom name"));
    }
    Residue_sp foundResidue = findResidue(molecule,residueIdentifier);
    if ( foundResidue.notnilp() )
    {
      if ( foundResidue->hasAtomWithName(atomName)) {
          LOG(BF("Found atom with name: %s") % _rep_(atomName)  );
	    Vector3 pos = foundResidue->atomWithName(atomName)->getPosition();
	    geom::OVector3_sp v = geom::OVector3_O::create();
	    v->setAll3(pos.getX(),pos.getY(),pos.getZ());
	    return v;
	}
      SIMPLE_ERROR(BF("Residue does not contain atom named: %s") % _rep_(atomName));
    }
    SIMPLE_ERROR(BF("Molecule does not contain residue with identifier: %s") % _rep_(residueIdentifier));
}







#define	FuncError(x) SIMPLE_ERROR(x)



/*! monomer 'id 'groupName :aliases ['comment]
 * Store the definition of a new monomer
 */



#define ARGS_chem__monomer "(monomerId groupName &optional monomerAliases comment)"
#define DECL_chem__monomer ""
#define DOCS_chem__monomer "monomer"
CL_LAMBDA(monomerId groupName &optional monomerAliases comment);
CL_DEFUN core::T_sp chem__monomer(core::Symbol_sp monomerId, core::Symbol_sp groupName, core::List_sp monomerAliases, core::String_sp comment)
{
    OligomerPart_Monomer_sp newMon;
    newMon = OligomerPart_Monomer_O::create();
    newMon->_MonomerId = monomerId;
    newMon->_GroupName = groupName;
    newMon->_MonomerAliases = monomerAliases;
    newMon->_Comment = "";
    if ( comment.notnilp() ) newMon->_Comment = comment->get_std_string();
    return newMon;
}

/*! link 'mon1id 'coupling mon2Cmd
 * Store the definition of a new link
 */



#define ARGS_chem__link "(monomer1Id coupling monomer2)"
#define DECL_chem__link ""
#define DOCS_chem__link "link"
CL_DEFUN core::T_sp chem__link(core::Symbol_sp monomer1Id, core::Symbol_sp coupling, OligomerPart_Monomer_sp monomer2 )
{
    OligomerPart_Link_sp newLink;
    newLink = OligomerPart_Link_O::create();
    newLink->_Monomer1Id = monomer1Id;
    newLink->_Coupling = coupling;
    newLink->_Monomer2 = monomer2;
    return newLink;
}



/*
__BEGIN_DOC( candoScript.general.calculatePoint, subsection, calculatePoint)
\scriptCmdRet{calculatePoint}{ String::distanceAtom Real::distance String::angleAtom Real::angleDegrees String::dihedralAtom Real::dihedralDegrees Object::positionLookup  }

Build a coordinate based on a distance, angle and dihedral angle from other atoms. The \sa{positionLookup} object can be either an IterateFragmentCoordinateEntries or a Residue.

__END_DOC
 */



#define ARGS_chem__calculate_point "()"
#define DECL_chem__calculate_point ""
#define DOCS_chem__calculate_point "calculatePoint"
CL_DEFUN core::T_sp chem__calculate_point()
{
    IMPLEMENT_ME(); // handle new way
#if 0
#if !KINEMATICS
    core::KeyedArguments_sp kargs = core::KeyedArguments_O::createFromKeyedObjectCons(args,lisp);
    string distanceAtom = kargs->getStringAndRemoveOrDefault("distanceAtom","");
    string angleAtom = kargs->getStringAndRemoveOrDefault("angleAtom","");
    string dihedralAtom = kargs->getStringAndRemoveOrDefault("dihedralAtom","");
    double distance = kargs->getDoubleAndRemoveOrDefault("distance",0.0);
    double angleDegrees = kargs->getDoubleAndRemoveOrDefault("angleDegrees",0.0);
    double dihedralDegrees = kargs->getDoubleAndRemoveOrDefault("dihedralDegrees",0.0);
    core::T_sp positionLookup = kargs->getAndRemove("positionLookup");

    // Calculate the position according to the arguments
    //
    Vector3 pnt,posDistanceAtom,posAngleAtom,posDihedralAtom;
    // If the 
    if ( positionLookup.isA<IterateFragmentCoordinateEntries_O>() )
    {
	IterateFragmentCoordinateEntries_sp it = positionLookup.as<IterateFragmentCoordinateEntries_O>();
	posDistanceAtom = it->positionOfAtomWithName(distanceAtom);
	posAngleAtom = it->positionOfAtomWithName(angleAtom);
	posDihedralAtom = it->positionOfAtomWithName(dihedralAtom);
    } else if ( positionLookup.isA<Residue_O>() )
    {
	Residue_sp res = positionLookup.as<Residue_O>();
	posDistanceAtom = res->positionOfAtomWithName(distanceAtom);
	posAngleAtom = res->positionOfAtomWithName(angleAtom);
	posDihedralAtom = res->positionOfAtomWithName(dihedralAtom);
    } else
    {
	SIMPLE_ERROR(BF("Unacceptable positionLookup class"));
    }
    LOG(BF("About to calculate point with values") );
    LOG(BF(" distance(%lf)  distancePoint(%s)") % distance % posDistanceAtom.asString().c_str()  );
    LOG(BF(" angleDegrees(%lf)  anglePoint(%s)") % angleDegrees % posAngleAtom.asString().c_str()  );
    LOG(BF(" dihedralDegrees(%lf)  dihedralPoint(%s)") % dihedralDegrees % posDihedralAtom.asString().c_str()  );
    pnt = buildUsingBondAngleDihedral(distance,posDistanceAtom,
				      angleDegrees*0.0174533,posAngleAtom,
				      dihedralDegrees*0.0174533,posDihedralAtom);
    LOG(BF("Generated point: %s") % pnt.asString().c_str()  );
    return geom::OVector3_O::createFromVector3(lisp,pnt);
#endif
#endif
}






CL_DEFUN core::T_sp chem__oligomer(core::List_sp parts)
{
    Oligomer_sp olig = Oligomer_O::create();
    core::HashTableEq_sp monomerMap = core::HashTableEq_O::create_default();
    CandoDatabase_sp bdb = getCandoDatabase();
    for ( auto p : parts ) {
	OligomerPart_Base_sp oligPart = p->car<OligomerPart_Base_O>();
	Monomer_sp mon = oligPart->createMonomer(bdb);
	olig->addMonomer(mon);
	monomerMap->setf_gethash(mon->getId(), mon);
	if ( oligPart.isA<OligomerPart_Link_O>() )
	{
          OligomerPart_Link_sp link = gc::As<OligomerPart_Link_sp>(oligPart);
	    core::Symbol_sp	mon1Id = link->_Monomer1Id;
	    core::Symbol_sp	mon2Id = link->_Monomer2->_MonomerId;
	    if ( !monomerMap->contains(mon1Id) ) SIMPLE_ERROR(BF("Unknown monomer id: %s")%_rep_(mon1Id));
	    if ( !monomerMap->contains(mon2Id) ) SIMPLE_ERROR(BF("Unknown monomer id: %s")%_rep_(mon2Id));
	    Monomer_sp mon1 = monomerMap->gethash(mon1Id).as<Monomer_O>();
	    Monomer_sp mon2 = monomerMap->gethash(mon2Id).as<Monomer_O>();
	    olig->couple(mon1,link->_Coupling,mon2);
	}
    }
    return olig;
}


extern "C" int gaff_debug;
CL_DEFUN void chem__set_gaff_debug(bool on) {
  gaff_debug = on;
}




void	setupCandoPrimitives(core::Lisp_sp env)
{

	// for Mbb package
//  Chem_temp_Defun(make_coordinate_array_from_atom_list);
//  Chem_temp_Defun(save_archive_with_auto_set_cando_database); // ,&prim_saveArchiveWithAutoSetCandoDatabase,env);
//  Chem_temp_Defun(standard_database); // , &prim_standardDatabase,env);
//  Chem_temp_Defun(database); // , &prim_database,ARGS_prim_database,env);
//    Defun(describeDatabase); // , &prim_describeDatabase,env);
//  Chem_temp_Defun(bundle_database_path); // , &prim_bundleDatabasePath,ARGS_prim_bundleDatabasePath,env);
//    defNoWrapPackage(ChemPkg,"dbval", &prim_dbval,env);
//  Chem_temp_Defun(monomer); // , &prim_monomer,env);
//  Chem_temp_Defun(link); // , &prim_link,env);
//  Chem_temp_Defun(load_mol2); // , &prim_loadMol2,env);
//  Chem_temp_Defun(save_mol2); // , &prim_saveMol2,env);
//  Chem_temp_Defun(atom_pos); // , &prim_atomPos,env);
//  Chem_temp_Defun(find_residue); // , &prim_findResidue,env);
//    Defun(setDatabase); // , &prim_setDatabase,env);
//    defNoWrapPackage(ChemPkg,"contextGrep", &prim_contextGrep,env);
//  Chem_temp_Defun(create_virtual_atom); // , &prim_createVirtualAtom,env);
//  Chem_temp_Defun(calculate_point); // , &prim_calculatePoint,env);
//    Defun(RequiredPlug); // , &prim_RequiredPlug,ARGS_prim_RequiredPlug,env);
};




void	setupPythonPrimitives(core::Lisp_sp env)
{

    IMPLEMENT_MEF("Add python primitives from candoScript.cc");
#if 0
	// for Mbb package
    boost::python::def_raw(ChemPkg,"saveArchiveWithAutoSetCandoDatabase",&prim_saveArchiveWithAutoSetCandoDatabase,"","",env);
    boost::python::def_raw(ChemPkg,"standardDatabase", &prim_standardDatabase,"","",env);
    boost::python::def_raw(ChemPkg,"database", &prim_database,ARGS_prim_database,"",env);
    boost::python::def_raw(ChemPkg,"describeDatabase", &prim_describeDatabase,"","",env);
    boost::python::def_raw(ChemPkg,"bundleDatabasePath", &prim_bundleDatabasePath,ARGS_prim_bundleDatabasePath,"",env);
//    boost::python::def_raw(ChemPkg,"dbval", &prim_dbval,"","",env);
    boost::python::def_raw(ChemPkg,"monomer", &prim_monomer,"","",env);
    boost::python::def_raw(ChemPkg,"link", &prim_link,"","",env);
    boost::python::def_raw(ChemPkg,"loadMol2", &prim_loadMol2,"","",env);
    boost::python::def_raw(ChemPkg,"saveMol2", &prim_saveMol2,"","",env);
    boost::python::def_raw(ChemPkg,"atomPos", &prim_atomPos,"","",env);
    boost::python::def_raw(ChemPkg,"findResidue", &prim_findResidue,"","",env);
    boost::python::def_raw(ChemPkg,"setDatabase", &prim_setDatabase,"","",env);
//    boost::python::def_raw(ChemPkg,"contextGrep", &prim_contextGrep,"","",env);
    boost::python::def_raw(ChemPkg,"createVirtualAtom", &prim_createVirtualAtom,"","",env);
    boost::python::def_raw(ChemPkg,"calculatePoint", &prim_calculatePoint,"","",env);
    boost::python::def_raw(ChemPkg,"RequiredPlug", &prim_RequiredPlug,ARGS_prim_RequiredPlug,"",env);
#endif
};








};
