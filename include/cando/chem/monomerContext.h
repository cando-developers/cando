#ifdef USE_TOPOLOGY
/*
    File: monomerContext.h
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


#ifndef MonomerContext_H  //[
#define MonomerContext_H
#include <clasp/core/common.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/symbolMap.h>
//#include	"conformation.h"
#include <cando/chem/atom.h>
#include <cando/chem/monomer.fwd.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/chem/residue.h>
//#include <cando/chem/candoDatabase.h>
#include <clasp/core/smallMap.h>
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
class MonomerContext_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,MonomerContext_O,"MonomerContext",core::CxxObject_O);
public:


public:
  typedef	gc::SmallMap<core::Symbol_sp,EntityNameSetBase_sp>	NeighborMap;
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
                   core::List_sp list,
                   MCStringStack& name);
public:

  CL_LISPIFY_NAME("setFocus");
  CL_DEFMETHOD 	void	setFocus(EntityNameSetBase_sp n )
  {
    this->_Self = n;
  };

  void	addNeighbor(core::Symbol_sp plugName, EntityNameSetBase_sp n);
  bool	hasNeighborWithCouplingName(core::Symbol_sp plugName ) { return this->_Neighbors.contains(plugName);};
  uint	numberOfNeighbors() { return this->_Neighbors.size(); };
  EntityNameSetBase_sp getNeighbor(core::Symbol_sp plugName) const;
  adapt::SymbolSet_sp getPlugNames() const;

  CL_LISPIFY_NAME("getFocus");
  CL_DEFMETHOD 	EntityNameSetBase_sp	getFocus() { return this->_Self; };

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
  core::HashTableEqual_sp	getAllSpecificKeys();
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
#endif
