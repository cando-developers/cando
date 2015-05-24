       
       

#ifndef	AtomIndexer_H //[
#define AtomIndexer_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/holder.h>
#include <cando/adapt/stringSet.fwd.h>

#include <cando/chem/chemPackage.h>

//#include "core/stringList.fwd.h"// atomIndexer.h wants StringList needs stringList.fwd.h

namespace chem {


SMART(AtomIndexer);
SMART(AtomIndexer);
class AtomIndexer_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,AtomIndexer_O,"AtomIndexer");

public:
//    void	archive(core::ArchiveP node);
private:
	vector<string>	_Names;
public:
	typedef	vector<string>::iterator	atomNameIterator;

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
	int	addAtomName(const string& name);
		/*! Set from a list of names separated by spaces */
	void	setFromAtomNames(const string& names);

		/*! Append a Cons of Text objects */
	void	appendConsOfTexts(core::List_sp list);

		/*! Set from a list of names in StringList*/
	void	setFromStringList(core::T_sp names);

	string	description() const;

	string	asString() const;

	int	numberOfAtomNames() { return this->_Names.size(); };
	bool	containsAtomName(const string& name);

	int	indexForAtomName(const string& name);



	DEFAULT_CTOR_DTOR(AtomIndexer_O);
};




SMART(MapOfMonomerNamesToAtomIndexers);
class MapOfMonomerNamesToAtomIndexers_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,MapOfMonomerNamesToAtomIndexers_O,"MapOfMonomerNamesToAtomIndexers");

public:
//    void	archive(core::ArchiveP node);
private:
    core::SymbolMap<AtomIndexer_O>	_AtomIndexers;
public:
    typedef	core::SymbolMap<AtomIndexer_O>::iterator	iterator;

public:
	iterator	begin() { return this->_AtomIndexers.begin(); };
	iterator	end() { return this->_AtomIndexers.end(); };

	core::SymbolSet_sp	getMonomerNamesAsSymbolSet() const;
	bool		recognizesMonomerName(core::Symbol_sp nm) const;

	AtomIndexer_sp	getAtomIndexerForMonomerName(core::Symbol_sp name) const;

	void set(core::Symbol_sp nm, AtomIndexer_O::smart_ptr ai);
	void eraseEntryForMonomer(core::Symbol_sp nm);



	DEFAULT_CTOR_DTOR(MapOfMonomerNamesToAtomIndexers_O);
};







};
TRANSLATE(chem::AtomIndexer_O);
TRANSLATE(chem::MapOfMonomerNamesToAtomIndexers_O);
#endif //]
