/*
    File: externalInterface.cc
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

#include <iostream>
#include <fstream>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "core/array.h"
#include "core/environment.h"
#include "structureComparer.h"
#include "superposableConformationCollection.h"
#include "externalInterface.h"
#include "conformationExplorer.h"
#include "core/wrappers.h"




namespace chem
{





//
// Constructor
//

    void	ExternalInterface_O::initialize()
    {
	this->Base::initialize();
	this->_Data = core::HashTableEq_O::create_default();
	this->_AtomNames = core::StringList_O::create();
	this->_WeakConformationExplorer = _Nil<ConformationExplorer_O>();
    }

//
// Destructor
//

#if INIT_TO_FACTORIES

#define ARGS_ExternalInterface_O_make "(file_name)"
#define DECL_ExternalInterface_O_make ""
#define DOCS_ExternalInterface_O_make "make ExternalInterface"
  ExternalInterface_sp ExternalInterface_O::make(ConformationExplorer_sp conformationExplorer, string stageName)
  {
    auto  me  = gctools::GC<ExternalInterface_O>::allocate_with_default_constructor();
    me->_WeakConformationExplorer = conformationExplorer;
    me->_StageName = stageName;
    return me;
  };

#else

    core::T_sp ExternalInterface_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	this->_WeakConformationExplorer = translate::from_object<ConformationExplorer_O>::convert(env->lookup(Pkg(),"conformationExplorer"));
	this->_StageName = translate::from_object<core::Str_O>::convert(env->lookup(Pkg(),"stageName"))->get();
	return _Nil<core::T_O>();
    }

#endif

#ifdef XML_ARCHIVE
    void	ExternalInterface_O::archiveBase(core::ArchiveP node)
    {
#ifdef	DEBUG_ON
	if ( node->saving() )
	{
	    ASSERTNOTNULLP(this->_WeakConformationExplorer,"WeakConformationExplorer is NULL");
	}
#endif
	node->archiveWeakPointer("WeakConformationExplorer",this->_WeakConformationExplorer);
	node->attribute("StageName",this->_StageName);
	node->archiveString("ExternalScript",this->_ExternalScript);
	node->attributeStringMap("ScriptSubstitutions",this->_ScriptSubstitutions);
	node->attribute("Data",this->_Data);
	node->attribute("Directory",this->_Directory);
	node->attribute("FileNamePrefix",this->_FileNamePrefix);
	node->attribute("AtomNames", this->_AtomNames);
	node->attributeStringMap("AtomNamesToAtoms",this->_AtomNamesToAtoms);
    }
#endif





    void	ExternalInterface_O::setConformationExplorer(ConformationExplorer_sp explorer)
    {
	ConformationCollection_O::atomIterator	ai;
	stringstream				atomName;
	this->_WeakConformationExplorer = explorer->sharedThis<ConformationExplorer_O>();
	ASSERTNOTNULL(this->_WeakConformationExplorer);
    }


/*! Assemble the filename for the gaussian file given the
 * conformation index and the file extension
 */
    boost_filesystem::path ExternalInterface_O::_assembleFileName( int entryIndex,
								    const string& stageName, const string& ext)
    {
	stringstream	fnStream;
	boost_filesystem::path	fileName;
	fnStream.str("");
	fnStream << this->_Directory << "/";
	fnStream << this->_FileNamePrefix
		 << "_" << this->externalInterfaceName() << "_"
		 << setw(4) << setfill('0') << entryIndex
		 << "_" << stageName << ext;
	LOG(BF("Assembled fileName: %s") % fnStream.str().c_str()  );
	fileName = fnStream.str();
	return fileName;
    }



/*! Assemble the filepath for the gaussian file given the
 * conformation index and the file extension
 */
    boost_filesystem::path ExternalInterface_O::_assembleFilePath( int entryIndex,
								    const string& stageName, const string& ext)
    {
	stringstream	fnStream;
	boost_filesystem::path	filePath;
	fnStream.str("");
	fnStream << this->_Directory << "/";
	fnStream << this->_FileNamePrefix
		 << "_" << this->externalInterfaceName() << "_"
		 << setw(4) << setfill('0') << entryIndex
		 << "_" << stageName << ext;
	LOG(BF("Assembled fileName: %s") % fnStream.str().c_str()  );
	filePath = fnStream.str();
	return filePath;
    }


    boost_filesystem::path ExternalInterface_O::assembleFileName(ConformationExplorerEntryStage_sp stage, const string& ext)
    {
	ConformationExplorerEntry_sp	entry;
	ConformationExplorer_sp		explorer;
	int				index;
	string				stageName;
	entry = stage->getConformationExplorerEntry();
	explorer = stage->getConformationExplorer();
	index = explorer->getEntryIndex(entry);
	stageName = stage->getStageName();
	LOG(BF("stageName(%s)") % stageName.c_str()  );
	return this->_assembleFileName(index,stageName,ext);
    }

    boost_filesystem::path ExternalInterface_O::assembleFilePath(ConformationExplorerEntryStage_sp stage, const string& ext)
    {
	ConformationExplorerEntry_sp	entry;
	ConformationExplorer_sp		explorer;
	int				index;
	string				stageName;
	entry = stage->getConformationExplorerEntry();
	explorer = stage->getConformationExplorer();
	index = explorer->getEntryIndex(entry);
	stageName = stage->getStageName();
	LOG(BF("stageName(%s)") % stageName.c_str()  );
	return this->_assembleFilePath(index,stageName,ext);
    }


    bool	ExternalInterface_O::parseSummaryFileForFinalConformationAndEnergy(ConformationExplorerEntryStage_sp entryStage)
    {
	ifstream	fin;
#define	BUFSIZE	2048
	char		buffer[BUFSIZE];
	string		line;
	vector<string>	parts;
	ConformationExplorer_sp	explorer;
	Atom_sp		atom;
	boost_filesystem::path	fileName;
	fileName = this->assembleFileName(entryStage,this->summaryFileExtension());
	LOG(BF("Final conformation and energy being parsed from a .sum (summary) file: %s") % fileName.string().c_str() );
	explorer = entryStage->getConformationExplorer();
	fin.exceptions(ios::goodbit);
	ASSERTP(boost_filesystem::exists(fileName),"The filename "+fileName.string()+" doesn't exist");
	ASSERTP(boost_filesystem::is_regular(fileName),"The filename "+fileName.string()+" isn't a regular file");
	fin.open(fileName.string().c_str(),ios::in);
	LOG(BF("Opened file") );
	uint	numCoordsRead = 0;
	while ( !fin.eof() )
	{
	    fin.getline(buffer,BUFSIZE);
	    line = buffer;
	    LOG(BF("Read line: %s") % line.c_str()  );
	    core::tokenize(line,parts," ");
	    if ( parts[0] == "coord" )
	    {
		string atomName = parts[1];
		Vector3 pos;
		double x = atof(parts[2].c_str());
		double y = atof(parts[3].c_str());
		double z = atof(parts[4].c_str());
		LOG(BF("Read coords: %s %lf %lf %lf") % atomName.c_str() % x % y % z  );
		pos.set(x,y,z);
		if ( this->_AtomNamesToAtoms.contains(atomName) )
		{
		    atom = this->_AtomNamesToAtoms.get(atomName);
		    entryStage->setCoordinateForAtom(atom,pos);
		} else
		{
		    SIMPLE_ERROR(BF("In summary file: "+fileName.string()+" did not recognize atom named: "+atomName));
		}
		numCoordsRead++;
	    } else if ( parts[0] == "hartrees" )
	    {
		double hartrees = atof(parts[1].c_str());
		LOG(BF("Read hartrees: %lf") % hartrees );
		entryStage->setEnergyKCal(hartrees*627.5095);
	    }
	}
	if ( numCoordsRead != explorer->numberOfAllAtoms() )
	{
	    SIMPLE_ERROR(BF("Not enough coordinates"));
	}
	return true;
    }







};
