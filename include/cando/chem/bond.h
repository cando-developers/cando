/*
    File: bond.h
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




#ifndef	BOND_H
#define	BOND_H

#include <clasp/core/common.h>
#include <cando/chem/bond.fwd.h>
#include <cando/chem/matter.h>

#include <cando/chem/constitutionAtoms.fwd.h>

#include <cando/chem/chemPackage.h>


namespace chem {






SMART(Bond);
class Bond_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,Bond_O,"Bond",core::CxxObject_O);
 public:
  bool fieldsp() const { return true;};
  void fields(core::Record_sp node);
 public:
  static void lisp_initGlobals(core::Lisp_sp lisp);

  friend	class Loop;
public:
  BondOrder	order;
  Atom_sp		_Atom1;
  Atom_sp		_Atom2;
    /*! Store Symbol keyed properties of bond in alist */
  core::List_sp	_Properties;
 public:
  void initialize();


 public:
	/*! Copy the bond but don't redirect atoms yet */
  Bond_sp copyDontRedirectAtoms();

	/*! Redirect atoms to their copy atoms */
  void addYourselfToCopiedAtoms();

 public:
  static Bond_sp create(Atom_sp a1, Atom_sp a2, BondOrder o= singleBond);
  static void canonicalizeBondOrder(Atom_sp& a1, Atom_sp& a2, BondOrder& order);
  static bool singleBondP(BondOrder order) {
    return (order == singleBond ||
            order == singleDashBegin ||
            order == singleDashEnd ||
            order == singleWedgeBegin ||
            order == singleWedgeEnd );
  };
 public:
	/*! Return true if atom is Atom1, false if its Atom2 and throw an exception if its neither */
  bool isAtom1(Atom_sp atom) const;

  virtual bool equalp(core::T_sp obj) const;

CL_LISPIFY_NAME("getAtom1");
CL_DEFMETHOD   Atom_sp getAtom1() const { return this->_Atom1; };
CL_LISPIFY_NAME("getAtom2");
CL_DEFMETHOD   Atom_sp getAtom2() const { return this->_Atom2; };

CL_LISPIFY_NAME("getOrder");
CL_DEFMETHOD   BondOrder getOrder()	{ return this->order; };
CL_DEFMETHOD int getOrderAsInt()	{ return this->order; };
  string	getOrderAsString();
  void	setOrder(BondOrder o)	{ this->order=o; };
CL_LISPIFY_NAME("setOrder");
CL_DEFMETHOD   void	setOrderFromInt(int o)	{ this->order=(BondOrder)(o); };
  Atom_sp	getOtherAtom(Atom_sp atom ) const;
  bool	isInterResidueBond(core::HashTable_sp atomToResidue);
  bool	isIntraResidueBond(core::HashTable_sp atomToResidue ) {return !this->isInterResidueBond(atomToResidue);};
  void	redirectToAtomCopies();

	/*! Dump all of the properties associated with this matter */
  string	propertiesAsString() const;
	    
		/*! Remove the property from this Matters property list
		 */
  void	clearProperty(core::Symbol_sp propertySymbol);
		/*! Define/set the value of the property.
		 * If it already exists it is overwritten.
		 */
  void	setProperty(core::Symbol_sp propertySymbol, core::T_sp value);
		/*! Define/set the value of the property to TRUE.
		 * If it already exists it is overwritten.
		 */
  core::T_sp getProperty(core::Symbol_sp propertySymbol, core::T_sp defVal );
		/*! Return true if the property exists.
		 */
  bool	hasProperty(core::Symbol_sp propertySymbol );


  bool isBondBetween(Atom_sp a, Atom_sp b) const;

  BondOrder maybeFlipOrder(Atom_sp from);


	/*! Return a ConstitutionBond_sp object using the ConstitutionBond_O::atomIdMap 
	  to map Atoms to ConstitutionAtomIndex0N values */
  ConstitutionBond_sp asConstitutionBond(Atom_sp from, const MapAtomsToConstitutionAtomIndex0N& atomMap);

		//! For bonds that are read from XML we instruct them to join the atoms that they are supposed to join
  void		joinYourAtoms();


  string description() const;
  string describeOther(Atom_sp from) const;

  bool		invalid(Atom_sp a);
  void		failIfInvalid(Atom_sp a);

		/*! Bonds that are archived from a stream need to tell their
		 * atoms that they are bonded to each other.
		 */
  void	imposeYourself();

  Bond_O( const Bond_O& bb );

 Bond_O() : _Properties(_Nil<core::T_O>()) {};
  virtual ~Bond_O() {};
};

extern	core::NullTerminatedEnumAssociation bondOrderKeys[];

string		bondOrderToString(BondOrder bo);
char            bondOrderToChar(BondOrder bo);
BondOrder	stringToBondOrder(string bos);

SMART(BondList);
class BondList_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,BondList_O,"BondList",core::CxxObject_O);
 public:
  bool fieldsp() const { return true;};
  void fields(core::Record_sp node);
 public:
  void initialize();
 public:
 private:
  gctools::Vec0_uncopyable<Bond_sp>	_Bonds;
//	gctools::Vec0<Bond_sp>	_Bonds;
 public:

  BondList_O(const BondList_O& );


 public:
//	typedef gctools::Vec0<Bond_sp>::iterator	iterator;
  typedef gctools::Vec0<Bond_sp>::iterator	iterator;
  typedef gctools::Vec0<Bond_sp>::const_iterator	const_iterator;

  iterator begin() { return this->_Bonds.begin(); };
  iterator end() { return this->_Bonds.end(); };
  Bond_sp		atIndex(int idx) { return this->_Bonds[idx];};
  void		addBond(Bond_sp b);
  void		append(Bond_sp b) { this->_Bonds.push_back(b);};
  void		removeBond(Bond_sp b);
  void		removeBondBetween(Atom_sp a,Atom_sp b);
  int		size() { return this->_Bonds.size(); };
  BondList_sp	deepishCopy() const;

  string	description() const;
  string describeOthers(Atom_sp from) const;

		//! Impose all of your bonds on their atoms (serialization)
  void	imposeYourself();


  DEFAULT_CTOR_DTOR(BondList_O);
};


};

TRANSLATE(chem::Bond_O);
TRANSLATE(chem::BondList_O);


namespace kw {
extern core::Symbol_sp& _sym_singleBond;
extern core::Symbol_sp& _sym_doubleBond;
extern core::Symbol_sp& _sym_tripleBond;
extern core::Symbol_sp& _sym_aromaticBond;
extern core::Symbol_sp& _sym_hydrogenBond;
extern core::Symbol_sp& _sym_virtualBond;
extern core::Symbol_sp& _sym_dashedSingleBond;
extern core::Symbol_sp& _sym_dashedDoubleBond;
extern core::Symbol_sp& _sym_singleWedgeBegin;
extern core::Symbol_sp& _sym_singleDashBegin;
extern core::Symbol_sp& _sym_singleWedgeEnd;
extern core::Symbol_sp& _sym_singleDashEnd;
};








#endif
