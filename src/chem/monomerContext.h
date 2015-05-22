       
       
//
// (C) 2004 Christian E. Schafmeister
//


#ifndef MonomerContext_H  //[
#define MonomerContext_H
#include <clasp/core/common.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/symbolSet.fwd.h>
#include <cando/chem/bond.h>
#include <cando/candoBase/vector3.h>
//#include	"conformation.h"
#include <cando/chem/atom.h>
#include <clasp/core/stringSet.fwd.h>
#include <cando/chem/residue.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/entityNameSet.h>


#include <cando/chem/chemPackage.h>


namespace chem {


SMART(SpecificContext);
SMART(SpecificContextSet);



#define	STRINGSTACKMAXPARTS	100
class MCStringStack {
private:
	string	parts[STRINGSTACKMAXPARTS];
	int	top;
public:
	void	clear() {this->top=0;};
	void	push(const string& s) {this->parts[this->top++] =s;HARD_ASSERT(this->top<STRINGSTACKMAXPARTS);};
	void	pop() {this->top--;};
	string	all() {
		stringstream	ss;
		int		i;
		ss.str("");
		for (i=0;i<=this->top;i++){
		    ss<<this->parts[i];
		}
		return ss.str();
	};
};




SMART(MonomerContext);
class MonomerContext_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,MonomerContext_O,"MonomerContext");
public:
//    void archiveBase(core::ArchiveP node);


public:
	typedef	core::SymbolMap<EntityNameSetBase_O>	NeighborMap;
private: // archive
	EntityNameSetBase_sp		_Self;
		//! Map coupling names to monomers
	NeighborMap		 	_Neighbors;

private:
	void	expandOutsSpecificContexts(
	    NeighborMap::iterator mit,
				SpecificContextSet_sp list,
				SpecificContext_sp name );
	void	expandOuts(NeighborMap::iterator mit,
				core::StringSet_sp list,
				MCStringStack& name);
public:

	void	setFocus(EntityNameSetBase_sp n )
	{
			this->_Self = n;
	};

	void	addNeighbor(core::Symbol_sp plugName, EntityNameSetBase_sp n);
	bool	hasNeighborWithCouplingName(core::Symbol_sp plugName ) { return this->_Neighbors.contains(plugName);};
	uint	numberOfNeighbors() { return this->_Neighbors.size(); };
	EntityNameSetBase_sp getNeighbor(core::Symbol_sp plugName) const;
	core::SymbolSet_sp getPlugNames() const;

	EntityNameSetBase_sp	getFocus() { return this->_Self; };

	/*! Return true if this MonomerContext recognizes
		this monomer in its context
	*/

	bool		containsMonomerContext(MonomerContext_sp mc);
	bool		recognizesMonomerAndEnvironment(Monomer_sp mon);

	void		setFocusFromMonomer(Monomer_sp mon);
    void		setFocusFromConstitutionName(core::Symbol_sp nm);

	bool		hasInCoupling();
	core::Symbol_sp	getInCoupling();
    void		testConsistency( std::ostream& out );

		/*! Each context expands to a set of SpecificContexts that
		 * list the specific context that this Context 
		 * recognizes.  This is done by expanding all 
		 * group names to monomer names. 
		 */
	SpecificContextSet_sp getAllSpecificContexts();
		/*! Each context expands to a set of strings that
		 * list the specific context that this Context 
		 * recognizes.  This is done by expanding all 
		 * group names to monomer names. 
		 */
	core::SymbolSet_sp	getAllSpecificKeys();
		//! Return the first specific key
    core::Symbol_sp		getFirstSpecificKey();

//	void	defineSpecificContextFromMonomer(Monomer_sp mon);


		/*! return a string that describes the context in a manner
			similar to getAllSpecificKeys except groups are not
			expanded to monomers
		*/
    core::Symbol_sp	getKey();

	string description() const;


	MonomerContext_O( const MonomerContext_O& tdb );
//	MonomerContext_O( CandoDatabase_sp tdb, Monomer_sp mon );



	DEFAULT_CTOR_DTOR(MonomerContext_O);
};



// Constructor function
extern MonomerContext_sp create_MonomerContext1();
//extern MonomerContext_sp create_MonomerContext2(Monomer_sp mon);
inline MonomerContext_sp create_MonomerContext()
{
    return create_MonomerContext1();
}








};
TRANSLATE(chem::MonomerContext_O);
#endif //]


