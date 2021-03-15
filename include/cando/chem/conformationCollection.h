/*
    File: conformationCollection.h
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
       
       
/*
 */

#ifndef ConformationCollection_H
#define	ConformationCollection_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/nVector.h>
#include <cando/geom/coordinateArray.h>



namespace       chem {


SMART(ConformationExplorer);
DEFINE_O_SMART_POINTERS(ConformationCollectionEntry);
DEFINE_O_SMART_POINTERS(ConformationCollection);
SMART(SuperposeEngine);


#define	ConformationCollection_AddedEntry	0


SMART(ConformationCollectionEntry);
class ConformationCollectionEntry_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,ConformationCollectionEntry_O,"ConformationCollectionEntry",core::CxxObject_O);

public:
  ConformationCollection_sp	_WeakConformationCollection;
  geom::SimpleVectorCoordinate_sp		_AllCoordinates;
  core::HashTableEq_sp		_Data;

public:
  void setAllCoordinates(geom::SimpleVectorCoordinate_sp ac);
  geom::SimpleVectorCoordinate_sp getAllCoordinates() { return this->_AllCoordinates; }
public:

  CL_DEFMETHOD ConformationCollection_sp	getConformationCollection() {_OF(); ASSERTNOTNULL(this->_WeakConformationCollection);return this->_WeakConformationCollection;};
  void	setConformationCollection(ConformationCollection_sp s);
//	core::HashTableEq_sp getData() { return this->_Data; };

  void	translateAllCoordinates(const Vector3& offset);

  void	writeCoordinatesToMatter(Matter_sp matter);
  void	extractCoordinatesFromMatter(Matter_sp matter);

  virtual	string	coordinatesAsString()	{_OF();SUBCLASS_MUST_IMPLEMENT();};
  virtual	bool	parseStructureFromFileName(const string& fileName) {_OF();SUBCLASS_MUST_IMPLEMENT();};

  void	initialize();



  DEFAULT_CTOR_DTOR(ConformationCollectionEntry_O);
};


#if 0
inline	ConformationCollectionEntry_sp	create_ConformationCollectionEntry(core::Lisp_sp e,ConformationCollection_sp s)
{
    GC_ALLOCATE(ConformationCollectionEntry_O, ee ); // RP_Create<ConformationCollectionEntry_O>(e);
    ee->setConformationCollection(s);
    return ee;
};
#endif




class ConformationCollection_O;
    typedef	gctools::smart_ptr<ConformationCollection_O>	ConformationCollection_sp;
SMART(ConformationCollection);
class ConformationCollection_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,ConformationCollection_O,"ConformationCollection",core::CxxObject_O);

public:
	void initialize();
public:
	Matter_sp				_Matter;
    gctools::Vec0<ConformationCollectionEntry_sp>	_Entries;
    gctools::SmallOrderedSet<Atom_sp>			_AllAtoms;
    core::HashTableEq_sp			_Data;
protected:

public:
//    static ConformationCollection_sp open_ConformationCollection(const string& fn);

public:
    typedef	gctools::Vec0<ConformationCollectionEntry_sp>::iterator	entryIterator;
    typedef	gctools::SmallOrderedSet<Atom_sp>::iterator	atomIterator;

public:
	geom::SimpleVectorCoordinate_sp	_SimpleVectorCoordinate(Matter_sp matter);

public:


    core::HashTableEq_sp	getData() { return this->_Data;};

	void	saveAs(const string& fn);


		/* Set the Matter object that this ConformationCollection
		 * compares conformations of
		 */
	virtual void	setMatter(Matter_sp agg);
	Matter_sp	getMatter() { return this->_Matter;};

	atomIterator begin_AllAtoms()
		{ return this->_AllAtoms.begin();};
	atomIterator end_AllAtoms()
		{ return this->_AllAtoms.end();};

	entryIterator	begin_Entries()
		{ return this->_Entries.begin();};
	entryIterator	end_Entries()
		{ return this->_Entries.end();};

  	entryIterator	begin_Entries() const
		{ return this->_Entries.begin();};
	entryIterator	end_Entries() const
		{ return this->_Entries.end();};

		//! Return the index of the entry
  size_t	getEntryIndex(ConformationCollectionEntry_sp entry);


  core::List_sp entriesAsList() const;
  
	virtual void	clearEntries();
	ConformationCollectionEntry_sp getEntry(uint i);
	ConformationCollectionEntry_sp firstEntry();

	void	translateAllCoordinates(const Vector3& offset);

		/* Return true if the conformation is a new one and
		 * should be added to the ConformationCollection
		 */

	/*! Add the entry to the ConformationCollection, return status
	 * Subclass this to add the entry in a sorted order
	 * or to not add entries if they don't meet some criteria
	 */
    virtual size_t addEntry(ConformationCollectionEntry_sp entry);

  void	addConformationCollection(ConformationCollection_sp structureList);
  void	addConformationExplorerLatestConformations(ConformationExplorer_sp structureList);
  void	addConformationExplorerSelectedStageConformations(ConformationExplorer_sp structureList);
  
  CL_DEFMETHOD size_t numberOfEntries() { return this->_Entries.size();};
  CL_DEFMETHOD size_t numberOfAllAtoms() { return this->_AllAtoms.size();};

	bool canRender() { return true; }
#ifdef RENDER
	geom::Render_sp rendered(core::List_sp options);
#endif

		/*! Create a new entry but dont add it yet
		 * Subclasses of ConformationCollection can return
		 * subclasses of ConformationCollectionEntry with different
		 * behavior
		 */
virtual	ConformationCollectionEntry_sp	createEntry();






DEFAULT_CTOR_DTOR(ConformationCollection_O);
};

};


TRANSLATE(chem::ConformationCollectionEntry_O);
TRANSLATE(chem::ConformationCollection_O);
#endif
