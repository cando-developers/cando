       
       
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

extern	string	ConformationCollectionEntry_Xml;



SMART(ConformationExplorer);
DEFINE_O_SMART_POINTERS(ConformationCollectionEntry);
DEFINE_O_SMART_POINTERS(ConformationCollection);
SMART(SuperposeEngine);


#define	ConformationCollection_AddedEntry	0


SMART(ConformationCollectionEntry);
class ConformationCollectionEntry_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ConformationCollectionEntry_O,"ConformationCollectionEntry");

public:
//    void	archive(core::ArchiveP node);
private:
	ConformationCollection_wp	_WeakConformationCollection;
        geom::CoordinateArray_sp		_AllCoordinates;
	core::HashTableEq_sp		_Data;

public:
        void setAllCoordinates(geom::CoordinateArray_sp ac);
        geom::CoordinateArray_sp getAllCoordinates() { return this->_AllCoordinates; }
public:

	ConformationCollection_sp	getConformationCollection() {_OF(); ASSERTNOTNULL(this->_WeakConformationCollection);return this->_WeakConformationCollection;};
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
class ConformationCollection_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ConformationCollection_O,"ConformationCollection");

public:
	void initialize();
public:
//	void	archiveBase(core::ArchiveP node);
public:
//	void	archive(core::ArchiveP node);
protected:
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
	geom::CoordinateArray_sp	_extractCoordinateArray(Matter_sp matter);

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

		//! Return the index of the entry
	uint	getEntryIndex(ConformationCollectionEntry_sp entry);


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
    virtual uint addEntry(ConformationCollectionEntry_sp entry);


	void	addConformationCollection(ConformationCollection_sp structureList);
	void	addConformationExplorerLatestConformations(ConformationExplorer_sp structureList);
	void	addConformationExplorerSelectedStageConformations(ConformationExplorer_sp structureList);

	uint		numberOfEntries() { return this->_Entries.size();};
	uint		numberOfAllAtoms() { return this->_AllAtoms.size();};

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
