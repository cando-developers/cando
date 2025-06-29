/*
    File: matter.h
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


/*
 *	matter.h
 *
 *	Superclass for Aggregate/Molecule/Residue/Atom
 *
 */

#ifndef Matter_H
#define	Matter_H
#include <clasp/core/common.h>

#include <iostream>
#include <vector>
#include <string>
#include <cando/geom/matrix.h>
//#include "geom/objectList.fwd.h"
#include <cando/geom/boundingBox.fwd.h>
#ifdef RENDER
//#include "geom/render.fwd.h"
#endif

#include <cando/chem/matter.fwd.h>
#include <cando/chem/atomIdMap.fwd.h>
#include <cando/chem/atomId.h>
#include <cando/chem/elements.h>
#include <clasp/core/hashTable.fwd.h>
#include <cando/chem/restraint.fwd.h>
//#include	"geom/vector3.h"


#define	MATTER_CLASS		99	//	'c'
#define	AGGREGATE_CLASS		103	//	'g'
#define	MOLECULE_CLASS		109	//	'm'
#define	RESIDUE_CLASS		114	//	'r'
#define	ATOM_CLASS		97	//	'a'
#define	VIRTUAL_ATOM_CLASS	118	//	'v'
#define	BOND_CLASS		98	//	'b'


#include <cando/chem/chemPackage.h>


namespace chem {

#define	UNDEFINED_ID -1




SMART(Alias);

typedef enum	MatterIdsEnum { containerId = MATTER_CLASS,
    			aggregateId = AGGREGATE_CLASS,
			moleculeId = MOLECULE_CLASS,
			residueId = RESIDUE_CLASS,
			atomId = ATOM_CLASS,
			bondId = BOND_CLASS } MatterId;


class	Loop;

SMART(Bond);
SMART(Atom);
SMART(Residue);
SMART(Molecule);
SMART(Aggregate);
SMART(Matter);





SMART(Matter );
class Matter_O : public core::CxxObject_O
{
  friend class ConstitutionAtoms_O;
  LISP_CLASS(chem,ChemPkg,Matter_O,"Matter",core::CxxObject_O);
 public:
  static Matter_sp make(string name);
 public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
//	void	serialize(serialize::SNode snode);

  friend	class Loop;
public:
  typedef gctools::Vec0_uncopyable<Matter_sp>::iterator contentIterator;
  typedef gctools::Vec0_uncopyable<Matter_sp>::const_iterator const_contentIterator;

protected:
  int			_Id;
  MatterName		_Name;
  /*! KEEP THIS as a vector
      A lot depends on residues
      maintaining an identical
      order of atoms
      Through the database
      building process */
  gctools::Vec0_uncopyable<Matter_sp>		_Contents;
  /*! Store Symbol keyed properties of matter
		 */
  core::List_sp   	_Properties;
	/*! Maintain a list of restraints that span this Matter_O object */
  core::List_sp	        _RestraintList;
 public:
	/*! Adjust the size of the contents array */
  void resizeContents(int sz);
	/*! Put a child at a particular content index */
  void putMatter( int index, Matter_sp matter );

 private:
  core::List_sp accumulateRestraints(core::List_sp allRestraints) const;
 public:

  friend class Aggregate_O;
  friend class Molecule_O;
  friend class Residue_O;
  friend class Atom_O;
  friend	Matter_sp Matter();
  friend	Matter_sp Matter(int i);

//	bool	equals(Matter_sp c);

  virtual string __repr__() const;
  virtual char getMatterType() { return MATTER_CLASS; };

  CL_LISPIFY_NAME("getId");
  CL_DEFMETHOD   int	getId() const;
  CL_DEFMETHOD   void setId(int id) { this->_Id = id; };
	/*! Accumulate all of the restraints in this matter and its contents into a single RestraintVector */
  core::List_sp allRestraints() const;

//  void	setTempFileId(int i) {this->_TempFileId = i;};
//  int	getTempFileId() { return this->_TempFileId; };

 public:
	/*! Return a deep copy of this Matter */
  virtual Matter_sp copy(core::T_sp new_to_old) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
 protected:
	/*! Internally used to create a copy of this Matter but not redirect Atoms to their copies.
	 Use copy->redirectAtoms() to direct the new matter to its copied atoms 
	once all the atoms have been copied.*/
  virtual Matter_sp copyDontRedirectAtoms(core::T_sp new_to_old);
	/*! Redirect atoms to their copies */
  virtual void redirectAtoms();

	/*! Copy the restraints from another object but dont' redirect atoms */
  virtual void copyRestraintsDontRedirectAtoms(Matter_sp orig);
	/*! Redirect restraint atoms */
  virtual void redirectRestraintAtoms();

 public:
  core::T_sp extract_coordinates(core::T_sp coords) const;
  Matter_sp apply_coordinates(core::T_sp coords);
  Matter_sp apply_coordinates_to_heavy_atoms(core::T_sp coords);
  contentIterator begin_contents() { return this->_Contents.begin(); };
  contentIterator end_contents() { return this->_Contents.end(); };
  const_contentIterator begin_contents() const { return this->_Contents.begin(); };
  const_contentIterator end_contents() const { return this->_Contents.end(); };
		/*! Transfer the coordinates from an equivalent (equal) Matter.
		 * If the (other) isnt equal then throw an exception.
		 */
  virtual void transferCoordinates(Matter_sp other) {_OF();SUBCLASS_MUST_IMPLEMENT();};

		/*! Return the first atom that matches the alias
		 */
#ifdef USE_TOPOLOGY
  virtual Atom_sp	aliasAtomOrNil(Alias_sp alias);
		/*! Return the first residue that matches the alias
		 */
  virtual Residue_sp aliasResidueOrNil(Alias_sp alias);
#endif
		/*! Return the first atom that matches the alias
		 */
  Atom_sp	aliasAtom(Alias_sp alias);
		/*! Return the first residue that matches the alias
		 */
  Residue_sp aliasResidue(Alias_sp alias);

  virtual bool applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value);
  void applyProperty(core::Symbol_sp prop, core::T_sp value);
  void applyPropertyList(core::List_sp list);
  
  void	setAllAtomMasks(int m);
  void	setAtomMasksForAllHeavyAtoms(int m);

	/*! Dump all of the properties associated with this matter */
  string	propertiesAsString() const;

		/*! Define/set the value of the property to TRUE.
		 * If it already exists it is overwritten.
		 */
  void	setPropertyTrue(core::Symbol_sp propertySymbol);

  CL_LISPIFY_NAME("properties");
  CL_DEFMETHOD   core::List_sp getProperties() { return this->_Properties; };
  CL_DEFMETHOD void setProperties(core::List_sp properties) { this->_Properties = properties; };
		/*! Remove the property from this Matters property list
		 */
  void	clearProperty(core::Symbol_sp propertySymbol);
		/*! Define/set the value of the property.
		 * If it already exists it is overwritten.
		 */
  void	setProperty(core::Symbol_sp propertySymbol, core::T_sp value);
		/*! Return the value of the property.
		 * Throw an exception if the property isn't defined.
		 */
  core::T_sp getProperty(core::Symbol_sp propertySymbol ) const;
		/*! Return the value of the property or the
		 * default if it isn't defined.
		 */
  core::T_sp getPropertyOrDefault(core::Symbol_sp propertySymbol, core::T_sp defVal );
		/*! Return true if the property exists.
		 */
  bool	hasProperty(core::Symbol_sp propertySymbol );

  CL_LISPIFY_NAME("numberOfAtoms");
  CL_DEFMETHOD   virtual uint	numberOfAtoms() const {_OF(); SUBCLASS_MUST_IMPLEMENT();};

  void	eraseContents(); // Empty the contents vector, don't free the memory

  contentIterator eraseContent(contentIterator x) {return this->_Contents.erase(x);};

  virtual void	makeAllAtomNamesInEachResidueUnique();
  virtual size_t fillInImplicitHydrogens();
  virtual void	randomizeAtomPositions();
  virtual void	perturbAtomPositions(double dist);



  CL_LISPIFY_NAME("setName");
  CL_DEFMETHOD   void	setName(MatterName sName) { this->_Name = sName; };
  MatterName getName() const { return this->_Name; };
  CL_LISPIFY_NAME("getName");
  CL_DEFMETHOD   MatterName getName_notConst() { return this->_Name; };

  virtual Matter_mv	addMatter( Matter_sp child );
  void	addMatterRetainId( Matter_sp child );
  void	removeMatter( Matter_sp child );

  Atom_sp		firstAtomWithName(MatterName name);

  bool          hasContent(Matter_sp maybe_child);
  bool		hasContentWithName(MatterName sName);
  Matter_sp	contentWithName(MatterName sName);
  core::T_sp	contentWithNameOrNil(MatterName sName);
  int		contentIndexWithName(MatterName sName);

  size_t        nextId() const;
  Matter_sp	contentWithId( int lid );
  bool		hasContentWithId( int lid );
  int		contentIndexWithId( int lid );

  core::List_sp		contentsAsList();
  core::List_sp		contents() { return this->contentsAsList();};

  CL_LISPIFY_NAME("isAggregate");
  CL_DEFMETHOD   virtual bool	isAggregate() { return false;}
  CL_LISPIFY_NAME("isMolecule");
  CL_DEFMETHOD   virtual bool	isMolecule() { return false;}
  CL_LISPIFY_NAME("isResidue");
  CL_DEFMETHOD   virtual bool	isResidue() {return false;}
  CL_LISPIFY_NAME("isAtom");
  CL_DEFMETHOD   virtual bool 	isAtom() { return false;}


		/*!
		 * Connect all atoms that are within 2A of each other
		 * that aren't already bonded to each other
		 */
  void connectAllCloseAtoms();

  virtual void clearRestraints();
  
	/*! Add a restraint to the matter */
  virtual void addRestraint(Restraint_sp restraint);

		/*! Invert the structure and any chiral or dihedral restraints
		 */
  virtual void	invertStructureAndRestraints();

  Matter_sp	contentWithStorageId( int lid );
  bool		hasContentWithStorageId( int lid );

  core::Vector_sp	allAtoms();
  core::Vector_sp	allAtomsOfElement(Element element);
  core::List_sp 		allAtomsOfElementAsList(Element element);

  void calculateVirtualAtomPositions();


//	void	set_StorageId(int sid) { this->_StorageId = sid; };
//	int	get_StorageId() { return this->_StorageId; };

  int		contentIndex( Matter_sp c);
  CL_LISPIFY_NAME("setContentAt");
  CL_DEFMETHOD   void setContentAt( size_t i, Matter_sp c );
  CL_LISPIFY_NAME("contentAt");
  CL_DEFMETHOD   Matter_sp	contentAt( size_t i ) const;
  CL_LISPIFY_NAME("contentSize");
  CL_DEFMETHOD   size_t		contentSize( ) { return this->_Contents.size(); };
  size_t length() const { return this->_Contents.size(); };
  void	translateAllAtoms(const Vector3& v);

  void	setAtomAliasesForResiduesNamed(core::List_sp residueAliasAtoms, core::List_sp atomAliases );

  virtual core::HashTable_sp atomToResidueMap();

  virtual void	applyTransformToAtoms( const Matrix& m );
  void applyTransformToRestraints(const Matrix& m);
  
  virtual	string	subMatter() {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
  virtual	int	totalNetResidueCharge();

  virtual	string	description() const { stringstream ss; ss << "container("<<_rep_(this->_Name)<<")@"<<std::hex<<this<<std::dec; return ss.str();}
  Vector3		geometricCenter();
  core::T_mv momentOfInertiaTensor();
  core::T_mv momentOfGeometry();
		/*! Return a bounding box for the matter padded with (pad).
		 */
  geom::BoundingCuboid_sp	matterCalculateBoundingCuboid(double pad);

  void		reparent(Matter_sp c);

#ifdef RENDER
  geom::Render_sp rendered(core::List_sp kopts);
#endif

	/*! Return all atoms as a Cons of atoms. */
  core::List_sp allAtomsAsList(bool allowVirtualAtoms) const;

	/*! Return all bonds as Cons of Cons (each entry is a Bond object) */
  core::List_sp allBondsAsList(bool allowVirtualAtoms) const;

	/*! Return all angles as Cons of Cons (each entry is an Angle) */
  core::List_sp allAnglesAsList(bool allowVirtualAtoms) const;

	/*! Return all proper torsions as Cons of Cons (each entry is ProperTorsion) */
  core::List_sp allProperTorsionsAsList(bool allowVirtualAtoms) const;

	/*! Return all improper torsions as Cons of Cons (each entry is an ImproperTorsion ) */
  core::List_sp allImproperTorsionsAsList(bool allowVirtualAtoms) const;

	/*! Build a map of AtomIds to Atoms */
  virtual AtomIdMap_sp buildAtomIdMap() const;

	/*! Return the atom with the AtomId */
  virtual Atom_sp atomWithAtomId(const AtomId& atomId) const;

	/*! Invert the stereochemistry of my immediate restraints */
  void invertStereochemistryOfRestraints();

  Matter_O(const Matter_O& c );
 Matter_O() :
//_NextContentId(1),
//    _Id(1),
//    _TempFileId(0),
  _Name(nil<core::Symbol_O>()),
  _Properties(nil<core::T_O>()),
  _RestraintList(nil<core::T_O>()) {};
};

};
#endif
