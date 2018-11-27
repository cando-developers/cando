/*
    File: specificContext.h
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
#ifndef	SpecificContext_H //[
#define SpecificContext_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <clasp/core/common.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/adapt/symbolSet.fwd.h>
#include <cando/chem/chemPackage.h>


namespace chem {


SMART(CandoDatabase);

struct	KeyEntry
{
    core::Symbol_sp	_Coupling;
    core::Symbol_sp	_Monomer;
};

SMART(SpecificContextSet);

SMART(SpecificContext );
class SpecificContext_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,SpecificContext_O,"SpecificContext",core::CxxObject_O);
    static SpecificContext_sp make();
public: // virtual functions inherited from Object
	void	initialize();
//	string	__repr__() const;

private: // instance variables
	core::Symbol_sp		_SelfName;
    gctools::Vec0<KeyEntry>	_Neighbors;
public:	// Creation class functions

public:

	void clear();

	void setSelfName(core::Symbol_sp selfName);

	void pushNeighbor(core::Symbol_sp coupling, core::Symbol_sp neighbor);
	void setLastNeighborName(core::Symbol_sp neighbor);
	void popNeighbor();

	bool allMonomersInDatabase(CandoDatabase_sp bdb);

	string asString();
    core::Symbol_sp asSymbol();

	SpecificContext_sp copy();

	SpecificContext_O( const SpecificContext_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(SpecificContext_O);
};






SMART(SpecificContextSet );
class SpecificContextSet_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,SpecificContextSet_O,"SpecificContextSet",core::CxxObject_O);
#if INIT_TO_FACTORIES
  public:
    static SpecificContextSet_sp make();
#else
    DECLARE_INIT();
#endif
public: // virtual functions inherited from Object
	void	initialize();
//	string	__repr__() const;

private: // instance variables
  gc::SmallMap<core::Symbol_sp, SpecificContext_sp>	_KeySets;

public:
  typedef gc::SmallMap<core::Symbol_sp,SpecificContext_sp>::iterator	iterator;

public:	// Creation class functions

protected:
    void addWithKey(core::Symbol_sp key, SpecificContext_sp keycontext);
public:

	SpecificContextSet_O::iterator begin() { return this->_KeySets.begin(); };
	SpecificContextSet_O::iterator end() { return this->_KeySets.end(); };

	core::List_sp asCons();

CL_LISPIFY_NAME("size");
CL_DEFMETHOD 	uint size() { return this->_KeySets.size(); };

	void add(SpecificContext_sp key);
	bool contains(SpecificContext_sp key);
	void remove(SpecificContext_sp key);
	void merge(SpecificContextSet_sp set);
		/*! Remove every member of (this) that is in (subSet)
		 */
	void remove(SpecificContextSet_sp subSet);

	adapt::SymbolSet_sp asSymbolSetOfKeys();

	SpecificContextSet_O( const SpecificContextSet_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(SpecificContextSet_O);
};

};
TRANSLATE(chem::SpecificContext_O);
TRANSLATE(chem::SpecificContextSet_O);
#endif //]
