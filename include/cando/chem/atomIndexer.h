/*
    File: atomIndexer.h
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
       
       

#ifndef	AtomIndexer_H //[
#define AtomIndexer_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/holder.h>
#include <cando/adapt/symbolSet.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/adapt/stringList.fwd.h>
#include <clasp/core/smallMap.h>

#include <cando/chem/chemPackage.h>

//#include "core/stringList.fwd.h"// atomIndexer.h wants StringList needs stringList.fwd.h

namespace chem {


SMART(AtomIndexer);
SMART(AtomIndexer);
class AtomIndexer_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,AtomIndexer_O,"AtomIndexer",core::CxxObject_O);

public:
//    void	archive(core::ArchiveP node);
private:
    gc::Vec0<core::Symbol_sp>	_Names;
public:
    typedef gc::Vec0<core::Symbol_sp>::iterator	atomNameIterator;

public:

	atomNameIterator begin_AtomNames()
		{ return this->_Names.begin(); };
	atomNameIterator end_AtomNames()
		{ return this->_Names.end(); };

//	bool	getIndexAll() { return this->_IndexAll; };
//	void	setIndexAll(bool b);

		/*! Add the name to the _Names list
		 * make sure that it is unique.
		 */
	int	addAtomName(core::Symbol_sp name);
		/*! Set from a list of names separated by spaces */
	void	setFromAtomNames(const string& names);

		/*! Append a Cons of Text objects */
	void	appendConsOfTexts(core::List_sp list);

		/*! Set from a list of names in StringList*/
	void	setFromStringList(adapt::StringList_sp names);

	string	description() const;

	string	asString() const;

	int	numberOfAtomNames() { return this->_Names.size(); };
	bool	containsAtomName(core::Symbol_sp name);

	int	indexForAtomName(core::Symbol_sp name);



	DEFAULT_CTOR_DTOR(AtomIndexer_O);
};




SMART(MapOfMonomerNamesToAtomIndexers);
class MapOfMonomerNamesToAtomIndexers_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,MapOfMonomerNamesToAtomIndexers_O,"MapOfMonomerNamesToAtomIndexers",core::CxxObject_O);

public:
//    void	archive(core::ArchiveP node);
private:
    adapt::SymbolMap<AtomIndexer_O>	_AtomIndexers;
public:
    typedef	adapt::SymbolMap<AtomIndexer_O>::iterator	iterator;

public:
	iterator	begin() { return this->_AtomIndexers.begin(); };
	iterator	end() { return this->_AtomIndexers.end(); };

	adapt::SymbolSet_sp	getMonomerNamesAsSymbolSet() const;
	bool		recognizesMonomerName(core::Symbol_sp nm) const;

	AtomIndexer_sp	getAtomIndexerForMonomerName(core::Symbol_sp name) const;

	void set(core::Symbol_sp nm, AtomIndexer_O::smart_ptr_type ai);
	void eraseEntryForMonomer(core::Symbol_sp nm);



	DEFAULT_CTOR_DTOR(MapOfMonomerNamesToAtomIndexers_O);
};







};
TRANSLATE(chem::AtomIndexer_O);
TRANSLATE(chem::MapOfMonomerNamesToAtomIndexers_O);
#endif //]
