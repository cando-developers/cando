/*
    File: externalInterface.h
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
       
       
#ifndef	ExternalInterface_H //[
#define ExternalInterface_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "boost/filesystem.hpp"
#include "core/common.h"
#include "conformationCollection.h"
#include "core/stringList.h"


#include "chemPackage.h"


namespace chem {


SMART(ConformationExplorer);
SMART(ConformationExplorerEntry);
SMART(ConformationExplorerEntryStage);
SMART(ZMatrix);
SMART(SuperposableConformationCollection);





SMART(ExternalInterface);
class ExternalInterface_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ExternalInterface_O,"ExternalInterface");
#if INIT_TO_FACTORIES
 public:
    static ExternalInterface_sp make(ConformationExplorer_wp conformationExplorer, string stageName);
#else
    DECLARE_INIT();
#endif
public:
	void initialize();
public:
public:
//	void	archiveBase(core::ArchiveP node);
private:
	ConformationExplorer_wp	_WeakConformationExplorer;
	string			_StageName;
	string			_ExternalScript;
	core::SymbolMap<core::Str_O>		_ScriptSubstitutions;
    core::HashTableEq_sp	_Data;
	string			_Directory;
	string			_FileNamePrefix;
protected:
	core::StringList_sp		_AtomNames;
	core::SymbolMap<Atom_O>		_AtomNamesToAtoms;
public:
	typedef	core::StringList_O::iterator	atomNameIterator;
protected:
	stringstream		_ErrorMessages;
protected:
	boost_filesystem::path _assembleFileName(int entryIndex, const string& stageName, const string& ext );
	boost_filesystem::path _assembleFilePath(int entryIndex, const string& stageName, const string& ext );
public:
	typedef	core::SymbolMap<core::Str_sp>::iterator	scriptSubstitutionIterator;
public:


	atomNameIterator	begin_AtomNames() { return this->_AtomNames->begin(); };
	atomNameIterator	end_AtomNames() { return this->_AtomNames->end(); };

		/*! Return the external interface name
		 */
	virtual	string	externalInterfaceName() {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
		/*! The extension for files that store the results
		 * that will be read back into CANDO: default .out
		 */
	virtual string	resultFileExtension() { return ".out"; };
		/*! Summary files are a special results file
		 * that summarizes the results, they can be created
		 * by external programs by parsing the external package
		 * results file
		 */
	virtual string	summaryFileExtension() { return ".sum";};

	virtual void	setConformationExplorer(ConformationExplorer_sp explorer);

	virtual boost_filesystem::path assembleFileName(ConformationExplorerEntryStage_sp stage, const string& extension );
	virtual boost_filesystem::path assembleFilePath(ConformationExplorerEntryStage_sp stage, const string& extension );


#ifdef	USEBOOSTPYTHON
	string	boost_assembleFileName(ConformationExplorerEntryStage_sp stage, const string& extension ) { return this->assembleFileName(stage,extension).string();};
	string	boost_assembleFilePath(ConformationExplorerEntryStage_sp stage, const string& extension ) { return this->assembleFilePath(stage,extension).string();};
#endif

	void	setFileNamePrefix(const string& prefix) { this->_FileNamePrefix = prefix; };
	void	setDirectory(const string& prefix) { this->_Directory = prefix; };
	string getDirectory() { return this->_Directory; };

	void	setExternalScript(const string& script) { this->_ExternalScript = script; };
	string	getExternalScript() { return this->_ExternalScript; };

	void	clearScriptSubstitutions();
	void	setScriptSubstitution(const string& key, const string& val);
	string	getScriptSubstitution(const string& key);
//	scriptSubstitutionIterator begin_ScriptSubstitutions() {return this->_ScriptSubstitutions.begin(); };
//	scriptSubstitutionIterator end_ScriptSubstitutions() {return this->_ScriptSubstitutions.end(); };


    core::HashTableEq_sp	getData()	{ return this->_Data;};

	string	getFileNamePrefix() { return this->_FileNamePrefix; };

	virtual string	getExternalCoordinatesAsTextForEntryStage(ConformationExplorerEntryStage_sp stage) {_OF();SUBCLASS_MUST_IMPLEMENT();};

	virtual bool	parseExternalResultFileForFinalConformationAndEnergy(ConformationExplorerEntryStage_sp entryStage) {_OF();SUBCLASS_MUST_IMPLEMENT();};

		/*! Results from QM calculations can be put in a ".sum" file and that can be used to bring results back into CANDO
		 * Each line of a summary file starts with 'coord' or 'hartrees'
		 */
	virtual bool	parseSummaryFileForFinalConformationAndEnergy(ConformationExplorerEntryStage_sp entryStage );

#if 0
	void	writeExternalInputFilesOptimizeNeedsMinimizerAtomsOnly();

	void	writeExternalInputFileWithZMatrix(ZMatrix_sp zMatrix);

#endif
	string	getErrorMessages() { return this->_ErrorMessages.str(); };


	ExternalInterface_O( const ExternalInterface_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(ExternalInterface_O);
};



};
#endif //]
