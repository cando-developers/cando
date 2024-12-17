/*
    File: trainer.h
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

#ifndef MolTrainers_H
#define	MolTrainers_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "core/holder.h"
#include "core/environment.fwd.h"

//#include "jobHistory.h"

namespace       chem
{

    SMART(CandoDatabase);
    SMART(MonomerContext);
    SMART(ConformationExplorer);

    class TrainerHeader_O;
    typedef	gctools::smart_ptr<TrainerHeader_O>	TrainerHeader_sp;



    SMART(TrainerHeader);
    SMART(TrainerHeader );
    class TrainerHeader_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,TrainerHeader_O,"TrainerHeader");

    public:
//	void	archive(core::ArchiveP node);
    public:
	static string trainerFileName() { return "trainer.cxml"; };
    public:
	string			_TopDirectory;
	core::Path_sp		_Path;
	string			_State;
	string			_StateMessage;
	/*! Save the maximum number of conformational/stochastic
	 * search fails that this trainer has seen
	 * This can be used to checkpoint stochastic searches
	 * and then restart them later until
	 * _MostConformationalSearchFailsSeen is >=
	 * _MaxConformationalSearchFails
	 */
	uint			_MostConformationalSearchFailsSeen;
	uint			_MaxConformationalSearchFailsRequired;
	Bignum	_SecondsSpentConformationalSearching;
	core::Binder_sp		_Data;
	MonomerContext_sp	_Context;
	uint			_UniqueIndex;
	TrainerHistory_sp	_History;
	/*! If the user removes trainers from the seed database
	 * then archived Trainers can find themselves no longer
	 * within the seed database.  This flag keeps track of that.
	 * We can keep these trainers around and not use them in
	 * the CoordinatesDatabase but if the seed ever expands
	 * again to contain them then they will be picked up and
	 * used again.
	 */
	bool			_NotPartOfSeedDatabase;


    public:
	void initialize();
    public:

	void setPath(core::Path_sp path) { this->_Path = path;};
	core::Path_sp getPath() { return this->_Path;};
	bool getNotPartOfSeedDatabase() { return this->_NotPartOfSeedDatabase;};
	void setNotPartOfSeedDatabase(bool b) { this->_NotPartOfSeedDatabase=b;};
	void setTopDirectory(const string& str)
	{
	    this->_TopDirectory = str;
	}
	string getTopDirectory()
	{
	    return this->_TopDirectory;
	}

	string generatePathName() const;

	void setSecondsSpentConformationalSearching(const Bignum& seconds);
	Bignum getSecondsSpentConformationalSearching();

	void resetMostConformationalSearchFailsSeen();
	void updateMostConformationalSearchFailsSeen(uint i);
	uint getMostConformationalSearchFailsSeen();
	/*! If there aren't enough conformational searches yet
	 * then return true
	 */
	bool conformationalSearchesRequired();

	void setMaxConformationalSearchFailsRequired(uint i);
	uint getMaxConformationalSearchFailsRequired();

	string getTrainerUniquePathName()
	{
	    stringstream ss;
	    ss << (fmt::format("trainer{:08d}" , this->getUniqueIndex() );
	    return ss.str();
	}

	core::Path_sp constructPath(string const& fileName ) const;




	/*! Construct the path to a file in the directory for this trainer
	  It will have the form path/fileName.extension */
	core::Path_sp constructFilePath( string const& fileName, string const& extension);

	void setState(const string& str);
	void setState(const string& str, const string& msg);
	void setState1(const string& str) { this->setState(str); };
	void setState2(const string& str,const string& msg) { this->setState(str,msg); };
	void setStateAndMessage(const string& str,const string& msg) { this->setState(str,msg); };
	string getState()
	{
	    return this->_State;
	}

	void setStateMessage(const string& str)
	{
	    this->_StateMessage = str;
	}
	string getStateMessage()
	{
	    return this->_StateMessage;
	}


	core::T_sp setData( core::Symbol_sp key, core::T_sp val);
	core::T_sp getData( core::Symbol_sp key );
	bool hasData( core::Symbol_sp key );
	MonomerContext_sp getContext();
	void setContext( MonomerContext_sp co);
	string getContextKey();

	void setUniqueIndex(uint i);
	uint getUniqueIndex();
	TrainerHistory_sp getHistory();
	core::Binder_sp getDataBinder();
	DEFAULT_CTOR_DTOR(TrainerHeader_O);
    };




    SMART(Trainer);
    class Trainer_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,Trainer_O,"Trainer");

    public:
//	void	archive(core::ArchiveP node);
    private:
	TrainerHeader_sp		_Header;
	core::Binder_sp			_Data;
    public:
	void initialize();

    public:
//	void setCandoDatabase(CandoDatabase_sp db);
	void setHeader(TrainerHeader_sp head);
	TrainerHeader_sp getHeader();
	core::T_sp setData( core::Symbol_sp key, core::T_sp data );
	bool hasData( core::Symbol_sp key);
	core::T_sp getData(core::Symbol_sp key);
	core::Binder_sp getDataBinder();

	//! Write the trainer to its file
	void write();
#if 0
	KeyValueIterator begin_Data() { return this->_Data->begin_keyValue(); };
	KeyValueIterator end_Data() { return this->_Data->end_keyValue(); };
	bool	canRender() { return true;};
#ifdef RENDER
	candoBase::Render_sp rendered(core::Cons_sp options);
#endif
#endif

    public:	// special functions for common entries
	ConformationExplorer_sp getConformationExplorer();
	void setConformationExplorer(ConformationExplorer_sp str );

	DEFAULT_CTOR_DTOR(Trainer_O);
    };





    core::HashTable_sp readTrainerHeader( core::Path_sp const& path);
    core::T_sp readTrainerHeader_cando( core::Cons_sp args,  core::Environment_sp environ, core::LispPtr lisp);


    extern Trainer_sp readTrainer(core::Path_sp const& path, core::LispPtr lisp);
    extern core::T_sp readTrainer_cando(core::Cons_sp args, core::Environment_sp environ, core::LispPtr lisp);




    SMART(TrainerOrganizer);
    class TrainerOrganizer_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,TrainerOrganizer_O,"TrainerOrganizer");

    public:
//	void	archive(core::ArchiveP node);

    private:
//	CandoDatabaseReference_sp 	_CandoDatabaseRef;
	core::StringMap<TrainerHeader_O>	_Entries;
	core::Binder_sp	_Data;
    public:
	typedef	core::StringMap<TrainerHeader_O>::iterator	entryIterator;
    public:
	void initialize();

    public:
	uint numberOfTrainerHeaders();
//	void setCandoDatabase(CandoDatabase_sp hold);
	void addTrainerHeader(TrainerHeader_sp head);
	TrainerHeader_sp getTrainerHeader(const string& key);
	bool hasTrainerHeaderWithTrainerContext(MonomerContext_sp context);

	/*! If this organizer has a TrainerHeader that recognizes
	 * a subset of the requested context then return it
	 * Otherwise return nil.
	 */
	TrainerHeader_sp getTrainerHeaderThatRecognizesSubsetOfContext(MonomerContext_sp context);

	core::T_sp setData( core::Symbol_sp key, core::T_sp data );
	bool hasData( core::Symbol_sp key);
	core::T_sp getData(core::Symbol_sp key);
	entryIterator begin_Entries() { return this->_Entries.begin();};
	entryIterator end_Entries() { return this->_Entries.end();};

	/*!Return a Cons list of pairs consisting of
	 * key/TrainerHeader
	 */
	core::Cons_sp entriesAsCons();

#if 0
	KeyValueIterator begin_Data() { return this->_Data->begin_keyValue(); };
	KeyValueIterator end_Data() { return this->_Data->end_keyValue(); };
#endif

    public:
	void setTopDirectory(const string& topDir);
	string getTopDirectory();
	void setCandoDatabaseFileName(const string& str);
	string getCandoDatabaseFileName();

	uint	getNextUniqueIndex();

	TrainerHeader_sp	findEntryUsingContextKey(const string& key);
	DEFAULT_CTOR_DTOR(TrainerOrganizer_O);
    };


    extern	void		writeTrainerOrganizer(TrainerOrganizer_sp org, core::LispPtr lisp);
    extern	TrainerOrganizer_sp 	readTrainerOrganizer(core::LispPtr lisp);








};


#endif
