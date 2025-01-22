/*
    File: trainer.cc
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

#include "core/common.h"
#include "trainer.h"
#include "core/array.h"
//#include "core/archive.h"
//#include "core/xmlSaveArchive.h"
//#include "core/xmlLoadArchive.h"
#include "core/bignum.h"
#include "core/fileSystem.h"
#include "geom/objectDictionary.h"
#include "monomerContext.h"
#include "conformationExplorer.h"
#include "core/wrappers.h"

#define DEFAULT_MAX_CONFORMATIONAL_SEARCH_FAILS_REQUIRED 60

namespace chem {

uint TrainerHeader_O::getUniqueIndex()
{
    return this->_UniqueIndex;
}

TrainerHistory_sp TrainerHeader_O::getHistory()
{
    return this->_History;
}
core::Binder_sp TrainerHeader_O::getDataBinder()
{
    return this->_Data;
}


void TrainerHeader_O::resetMostConformationalSearchFailsSeen()
{
    this->_MostConformationalSearchFailsSeen = 0;
}

void TrainerHeader_O::updateMostConformationalSearchFailsSeen(uint i)
{
    if ( i > this->_MostConformationalSearchFailsSeen )
    {
	this->_MostConformationalSearchFailsSeen = i;
    }
}

uint TrainerHeader_O::getMostConformationalSearchFailsSeen()
{
    return this->_MostConformationalSearchFailsSeen;
}

bool TrainerHeader_O::conformationalSearchesRequired()
{
    return this->_MostConformationalSearchFailsSeen<this->_MaxConformationalSearchFailsRequired;
}


void TrainerHeader_O::setMaxConformationalSearchFailsRequired(uint i)
{
    this->_MaxConformationalSearchFailsRequired = i;
}

uint TrainerHeader_O::getMaxConformationalSearchFailsRequired()
{
    return this->_MaxConformationalSearchFailsRequired;
}


void TrainerHeader_O::setSecondsSpentConformationalSearching(const Bignum& i)
{
    this->_SecondsSpentConformationalSearching = i;
}

Bignum TrainerHeader_O::getSecondsSpentConformationalSearching()
{
    return this->_SecondsSpentConformationalSearching;
}





core::Path_sp TrainerHeader_O::constructPath( string const& fileName ) const
    {
	core::Path_sp path = core::Path_O::create(this->_Path->getPath());
	core::Path_O* p = path.get();
	(*p) /= fileName;
	return path;
    }

    core::T_sp TrainerHeader_O::setData( core::Symbol_sp key, core::T_sp val)
{
    this->_Data->extend(key,val);
    return val;
}

    core::T_sp TrainerHeader_O::getData( core::Symbol_sp key )
{
    return this->_Data->lookup(key);
}

    bool TrainerHeader_O::hasData( core::Symbol_sp key )
{
    return this->_Data->contains(key);
}

MonomerContext_sp TrainerHeader_O::getContext()
{
    return this->_Context;
}
void TrainerHeader_O::setContext( MonomerContext_sp co)
{
    this->_Context = co;
}

void	TrainerHeader_O::setState(const string& str)
{
    this->_State = str;
}

void TrainerHeader_O::setState(const string& str, const string& msg)
{
    this->setState(str);
    this->setStateMessage(msg);
}



void	TrainerHeader_O::setUniqueIndex(uint i)
{
    this->_UniqueIndex = i;
}

void	TrainerHeader_O::initialize()
{
    this->Base::initialize();
    this->_Data = core::Binder_O::create();
    this->_Context = _Nil<MonomerContext_O>();
    this->_History = TrainerHistory_O::create();
    this->_TopDirectory = "";
    this->_State = "";
    this->_StateMessage = "";
    this->_History->clear();
    this->_Context = _Nil<MonomerContext_O>();
    this->_MostConformationalSearchFailsSeen = 0;
    this->_MaxConformationalSearchFailsRequired = DEFAULT_MAX_CONFORMATIONAL_SEARCH_FAILS_REQUIRED;
    this->_SecondsSpentConformationalSearching = 0;
}



#ifdef XML_ARCHIVE
    void	TrainerHeader_O::archive(core::ArchiveP node)
{
    node->attribute("TopDirectory", this->_TopDirectory );
    node->attribute("Path", this->_Path );
    node->attribute("State", this->_State );
    node->archiveStringIfNotDefault("StateMessage", this->_StateMessage, "" );
    node->attributeIfNotDefault("NotPartOfSeedDatabase", this->_NotPartOfSeedDatabase, true);
    node->attribute("data",this->_Data);
    node->attribute("context",this->_Context);
    node->attribute("UniqueIndex", this->_UniqueIndex );
    node->attribute("history",this->_History);
    node->attributeIfNotDefault("MostConformationalSearchFailsSeen",this->_MostConformationalSearchFailsSeen,(uint)(0));
    node->attributeIfNotDefault("MaxConformationalSearchFailsRequired",this->_MaxConformationalSearchFailsRequired,(uint)(DEFAULT_MAX_CONFORMATIONAL_SEARCH_FAILS_REQUIRED));
    node->attribute("SecondsSpentConformationalSearching",this->_SecondsSpentConformationalSearching);
}
#endif

#if 0
string TrainerHeader_O::generateFileName() const
{
    LOG("Generating file name directory(%s)" , this->_TopDirectory.c_str() );
    LOG("Generating file name fileNamePrefix(%s)" , this->getTrainerFileNamePrefix().c_str()  );
    stringstream fileName;
    fileName << (fmt::format("{}/{}" , this->_TopDirectory , this->getTrainerFileNamePrefix() );
    LOG("Put it all together fileName(%s)" , fileName.str()  );
    return fileName.str();
}
#endif

string TrainerHeader_O::getContextKey()
{
    return this->_Context->getKey();
}





void	Trainer_O::initialize()
{
    this->Base::initialize();
    this->_Data = core::Binder_O::create();
};

#ifdef XML_ARCHIVE
    void	Trainer_O::archive(core::ArchiveP node)
{
    node->attribute("header",this->_Header);
    node->attribute("data",this->_Data);
}
#endif

void	Trainer_O::write()
{
    IMPLEMENT_ME();
#ifdef XML_ARCHIVE
    core::XmlSaveArchive_sp a = core::XmlSaveArchive_O::create();
    core::Path_sp path = this->getHeader()->constructPath(this->getHeader()->trainerFileName());
    a->put("job",this->sharedThis<Trainer_O>());
    a->saveAs(path->asString());
#endif
}

ConformationExplorer_sp Trainer_O::getConformationExplorer()
{
  return (this->getData(_lisp->internWithPackageName(ChemPkg,"_ConformationExplorer"))).as<ConformationExplorer_O>();
};
void Trainer_O::setConformationExplorer(ConformationExplorer_sp str )
{
    this->setData(_lisp->internWithPackageName(ChemPkg,"_ConformationExplorer"),str);
};

#if 0
geom::Render_sp	Trainer_O::rendered(core::Cons_sp opts)
{
    return this->_Data->rendered(opts);
}
#endif
void Trainer_O::setHeader(TrainerHeader_sp head)
{
    this->_Header = head;
}
TrainerHeader_sp Trainer_O::getHeader()
{
    return this->_Header;
}

    core::T_sp Trainer_O::setData( core::Symbol_sp key, core::T_sp data )
{
    this->_Data->extend(key,data);
    return data;
}
    bool Trainer_O::hasData( core::Symbol_sp key)
{
    return this->_Data->contains(key);
}
    core::T_sp Trainer_O::getData(core::Symbol_sp key)
{
    return this->_Data->lookup(key);
}

    core::Binder_sp Trainer_O::getDataBinder()
{
    return this->_Data;
}


core::HashTable_sp readTrainerHeader(core::Path_sp path, core::LispPtr lisp )
{
    char buffer[1024];
    string fileName = path->asString();
    ifstream fin(fileName.c_str());
    if ( fin.bad() )
    {
	SIMPLE_ERROR(("Could not open file: "+fileName));
    }
    vector<string> parts;
    uint lineCount = 0;
    core::Instance_sp trainerHeaderClass = lisp->classFromClassName(TrainerHeader_O::static_classSymbol());
    string search = trainerHeaderClass->getPackagedName();
    while ( !fin.eof() )
    {
	lineCount++;
	if ( lineCount > 100 )
	{
          SIMPLE_ERROR(("Could not find (%s) in the first 100 lines") , search);
	}
	fin.getline(buffer,1024);
	string line = buffer;
	string trimmed = core::trimWhiteSpace(line);
	if ( trimmed.size() < search.size() ) continue;
	if ( trimmed.substr(1,search.size()) == search )
	{
	    trimmed = trimmed.substr(0,trimmed.size()-1); // remove ">"
	    parts = core::split(trimmed);
	    break;
	}
    }
    fin.close();
    vector<string>::iterator it;
    it = parts.begin(); it++; // Throw away the Html Line
    core::HashTable_sp dict = core::HashTable_O::createEq();
    for ( ; it<parts.end(); it++ )
    {
	vector<string> onePair = core::split(*it,"=");
	dict->put(onePair[0],core::Str_O::create(onePair[1].substr(1,onePair[1].size()-2)));
    }
    return dict;
}




    
    
#define ARGS_af_readTrainerHeader "(path)"
#define DECL_af_readTrainerHeader ""
#define DOCS_af_readTrainerHeader "readTrainerHeader"
core::T_sp af_readTrainerHeader(core::Path_sp path)
    {
	return readTrainerHeader(path,_lisp);
    }




#ifdef XML_ARCHIVE
Trainer_sp readTrainer(core::Path_sp const& path, core::LispPtr lisp)
{
    Trainer_sp job;
    core::XmlLoadArchive_sp a;
    a = core::XmlLoadArchive_O::create();
    core::Path_sp trainerPath = path;
    trainerPath->path_append(TrainerHeader_O::trainerFileName());
    string fileName = trainerPath->asString();
    a->parse(fileName);
    job = (a->get("job")).as<Trainer_O>();
    return job;
}
#endif


#if 0
Trainer_sp readTrainerReplaceDatabase(const string& fileName, CandoDatabase_sp bdh)
{
    Trainer_sp job;
    TrainerHeader_sp jobHeader;
    job = readTrainer(fileName);
    jobHeader = job->getHeader();
    		//
		// When we change databases we will break the link
		// between the FocusResidue and its constitution
		// Swap it here
    Residue_sp focusResidue = (job->getData("FocusResidue")).as<Residue_O>();
    Constitution_sp focusResidueConstitution = focusResidue->getConstitution();
    string constitutionName = focusResidueConstitution->getName();
    chem::CandoDatabase_sp oldTrainerCandoDatabase = job->getCandoDatabase();
    LOG("Replacing AlchemistDatabse of job with new one" );
    oldTrainerCandoDatabase->giveYourDependantObjectsTo(bdh);
    Constitution_sp newFocusResidueConstitution = bdh->monomerConstitutionForName(constitutionName);
    LOG("Setting constitution for FocusResidue to: %s" , newFocusResidueConstitution->description().c_str()  );
    focusResidue->setConstitution(newFocusResidueConstitution);
    LOG("Setting constitution for FocusResidue: %s" , focusResidue->description().c_str()  );
    LOG("    constitution = %s" , newFocusResidueConstitution->description().c_str()  );
    job->setCandoDatabase(bdh);
    return job;
}
#endif



#if 0
    void writeTrainer(Trainer_sp job, const string& pathName, core::LispPtr lisp)
{
    core::XmlSaveArchive_sp a;
    job->getHeader()->setPathName(pathName);
    a = job->lisp()->create<core::XmlSaveArchive_O>();
    a->put("job",job);
    LOG("writeTrainer to file: %s" , fileName.c_str()  );
    a->saveAs(fileName);
}
#endif







void TrainerOrganizer_O::initialize()
{
    this->Base::initialize();
    this->_Data = core::Binder_O::create();
}


#ifdef XML_ARCHIVE
void	TrainerOrganizer_O::archive(core::ArchiveP node)
{
    node->attributeStringMap("entries",this->_Entries);
    node->attribute("data",this->_Data);
}
#endif

uint TrainerOrganizer_O::numberOfTrainerHeaders() 
{ 
    return this->_Entries.size(); 
};

    core::T_sp TrainerOrganizer_O::setData( core::Symbol_sp key, core::T_sp data )
{
    this->_Data->extend(key,data);
    return data;
}
bool TrainerOrganizer_O::hasData( core::Symbol_sp key)
{
    return this->_Data->contains(key);
}
    core::T_sp TrainerOrganizer_O::getData(core::Symbol_sp key)
{
    return this->_Data->lookup(key);
}


#if 0
    core::StringMap<TrainerHeader_O>::value_iterator TrainerOrganizer_O::begin_Entries() { return this->_Entries.begin_value(); };
    core::StringMap<TrainerHeader_O>::value_iterator TrainerOrganizer_O::end_Entries() { return this->_Entries.end_value(); };
#endif

void TrainerOrganizer_O::setTopDirectory(const string& topDir)
{
    core::Str_sp val;
    val = core::Str_O::create(topDir);
   this->_Data->extend(_lisp->internWithPackageName(ChemPkg,"_TopDirectory"), val );
}
string TrainerOrganizer_O::getTopDirectory()
{
  return (this->_Data->lookup(_lisp->internWithPackageName(ChemPkg,"_TopDirectory"))).as<core::Str_O>()->get();
}
void TrainerOrganizer_O::setCandoDatabaseFileName(const string& str)
{
    core::Str_sp val;
    val = core::Str_O::create(str);
    this->_Data->extend(_lisp->internWithPackageName(ChemPkg,"_CandoDatabaseFileName"), val );
}
string TrainerOrganizer_O::getCandoDatabaseFileName()
{
  return (this->_Data->lookup(ChemPkg,"_CandoDatabaseFileName")).as<core::Str_O>()->get();
}


void TrainerOrganizer_O::addTrainerHeader(TrainerHeader_sp head)
{
    this->_Entries.set(head->getContextKey(),head);
}
TrainerHeader_sp TrainerOrganizer_O::getTrainerHeader(const string& key)
{
    return this->_Entries.get(key);
}

bool	TrainerOrganizer_O::hasTrainerHeaderWithTrainerContext(MonomerContext_sp context )
{
string	key;
    key = context->getKey();
    return this->_Entries.contains(key);
}




TrainerHeader_sp	TrainerOrganizer_O::findEntryUsingContextKey(const string& key)
{
    ASSERT(this->_Entries.contains(key));
    return this->_Entries.get(key);
}



TrainerHeader_sp	TrainerOrganizer_O::getTrainerHeaderThatRecognizesSubsetOfContext(MonomerContext_sp context)
{
    for ( entryIterator ei=this->_Entries.begin(); 
    		ei!=this->_Entries.end(); ei++ )
    {
        MonomerContext_sp organizerContext = (ei->second)->getContext();
	if ( context->containsMonomerContext(organizerContext) )
	{
	    return (ei->second);
	}
    }
    return _Nil<TrainerHeader_O>();
}


uint	TrainerOrganizer_O::getNextUniqueIndex()
{
entryIterator	ei;
uint		highestUniqueIndex = 0;
    for ( ei=this->_Entries.begin(); ei!=this->_Entries.end(); ei++ )
    {
	highestUniqueIndex = max(highestUniqueIndex,(ei->second)->getUniqueIndex());
    }
    if ( highestUniqueIndex == UndefinedUnsignedInt ) return 0;
    return highestUniqueIndex+1;
}



core::Cons_sp	TrainerOrganizer_O::entriesAsList()
{
    core::StringMap<TrainerHeader_O>::iterator	it;
    core::Cons_sp entries, cur, tail;
    entries = _Nil<core::Cons_O>();
    for ( it=this->_Entries.begin(); it!=this->_Entries.end(); it++ )
    {
	core::Cons_sp one = core::Cons_O::createList(core::Str_O::create(it->first),it->second);
	cur = core::Cons_O::create(one);
	if ( entries.nilp() )
	{
	    entries = cur;
	    tail = cur;
	} else
	{
	    tail->setCdr(cur);
	    tail = cur;
	}
    }
    return entries;
}

#ifdef XML_ARCHIVE
    void	writeTrainerOrganizer(TrainerOrganizer_sp org,core::LispPtr lisp)
{
    core::XmlSaveArchive_sp	archive;
    archive = org->lisp()->create<core::XmlSaveArchive_O>();
    archive->put("organizer",org);
    archive->saveAs("_organizer.cxml");
}


TrainerOrganizer_sp readTrainerOrganizer(core::LispPtr lisp)
{
    core::XmlLoadArchive_sp	archive;
    TrainerOrganizer_sp org;
    archive = core::XmlLoadArchive_O::create();
    archive->parse("_organizer.cxml");
    org = (archive->get("organizer")).as<TrainerOrganizer_O>();
    return org;
}







#define ARGS_af_readTrainer "(path)"
#define DECL_af_readTrainer ""
#define DOCS_af_readTrainer "readTrainer"
core::T_sp af_readTrainer(core::Path_sp path)
{
    Trainer_sp train = readTrainer(path,_lisp);
    return train;
}
#endif




















};
