#ifdef USE_TOPOLOGY

/*
    File: monomerPack.h
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
       
       
#ifndef	MonomerPack_H
#define MonomerPack_H




#include <clasp/core/foundation.h>
#include <cando/chem/entityNameSet.h>
#include <cando/chem/atomIndexer.h>


#include <cando/chem/chemPackage.h>

#include <clasp/core/cons.fwd.h>// monomerPack.h wants Cons needs cons.fwd.h
#include <cando/adapt/symbolList.fwd.h>// entityNameSet.h wants SymbolList needs symbolList.fwd.h

namespace chem {

SMART(MapOfMonomerNamesToAtomIndexers);
SMART(Cons);
SMART(SymbolList);

SMART(MonomerPack);
class MonomerPack_O : public EntityNameSet_O
{
    LISP_CLASS(chem,ChemPkg,MonomerPack_O,"MonomerPack",EntityNameSet_O);

public:
	void initialize();
public:
		//! Store the unique IDs of the Interesting atoms for each monomer
        core::List_sp				_InterestingAtomAliases;
	MapOfMonomerNamesToAtomIndexers_sp	_AtomIndexers;

private:
	void _checkAtomIndexers();

public:

//	void	setInterestingAtomAliasesFromString(const string& nm);
	void	setInterestingAtomAliasesFromSymbolList(adapt::SymbolList_sp nm);
	string	getInterestingAtomAliasesAsString();

        virtual core::List_sp getInterestingAtomAliases() { return this->_InterestingAtomAliases; };

	void	setMonomerNameOrPdb(core::Symbol_sp nm);
	void	addMonomerName(core::Symbol_sp nm);
	void	removeMonomerName(core::Symbol_sp nm);

		/*! Build the contents from information in a Cons */
	void	defineContentsFromCons(core::List_sp atomAliases, core::List_sp parts);

	virtual string	getMonomerNameWithAtoms(core::Symbol_sp nm);

	void extendAliases( core::List_sp atomAliases, core::List_sp parts);

    virtual	bool	supportsInterestingAtomAliases() { return true;};
    virtual void setInterestingAtomNamesForMonomerName(core::Symbol_sp monomerName, const string& atomIndexerNames);
  virtual void setInterestingAtomNamesForMonomerNameStringList(core::Symbol_sp monomerName, core::List_sp atomIndexerNames);
    virtual void setInterestingAtomNamesForMonomerNameFromCons(core::Symbol_sp, core::List_sp atomIndexerNames);
    virtual string getInterestingAtomNamesForMonomerName(core::Symbol_sp name);
    virtual bool hasInterestingAtomAlias(Alias_sp alias);
    virtual int getInterestingAtomAliasIndex(Alias_sp alias);
    virtual AtomIndexer_sp getAtomIndexerForMonomerName(core::Symbol_sp monomerName);

	MonomerPack_O( const MonomerPack_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(MonomerPack_O);
};



};


TRANSLATE(chem::MonomerPack_O);
#endif
#endif
