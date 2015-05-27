#ifndef	SpecificContext_H //[
#define SpecificContext_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <clasp/core/common.h>
#include <clasp/core/holder.h>
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
class SpecificContext_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,SpecificContext_O,"SpecificContext");
    static SpecificContext_sp make();
public: // virtual functions inherited from Object
	void	initialize();
//	void	archive(core::ArchiveP node);
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
class SpecificContextSet_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,SpecificContextSet_O,"SpecificContextSet");
#if INIT_TO_FACTORIES
  public:
    static SpecificContextSet_sp make();
#else
    DECLARE_INIT();
#endif
public: // virtual functions inherited from Object
	void	initialize();
//	void	archive(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
	adapt::SymbolMap<SpecificContext_O>	_KeySets;

public:
	typedef adapt::SymbolMap<SpecificContext_O>::iterator	iterator;

public:	// Creation class functions

protected:
    void addWithKey(core::Symbol_sp key, SpecificContext_sp keycontext);
public:

	SpecificContextSet_O::iterator begin() { return this->_KeySets.begin(); };
	SpecificContextSet_O::iterator end() { return this->_KeySets.end(); };

	core::List_sp asCons();

	uint size() { return this->_KeySets.size(); };

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
