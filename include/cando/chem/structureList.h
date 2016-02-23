/*
    File: structureList.h
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

#ifndef StructureList_H
#define	StructureList_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
//#include	"conformation.h"
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/nVector.h>
#include <cando/geom/coordinateArray.h>



namespace       chem {

extern	string	StructureListEntry_Xml;

DEFINE_O_SMART_POINTERS(Structure_Old_ListEntry);
DEFINE_O_SMART_POINTERS(Structure_Old_List);


#define	StructureList_AddedEntry	0


SMART(StructureListEntry);
class Structure_Old_ListEntry_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,Structure_Old_ListEntry_O,"StructureListEntry",core::CxxObject_O);

public:
//    void	archive(core::ArchiveP node);
	void	initialize();
public:
static Structure_Old_ListEntry_sp	create(Structure_Old_List_sp s);
public:
	Structure_Old_List_wp		_WeakStructureList;
        geom::CoordinateArray_sp	_AllCoordinates;
        geom::CoordinateArray_sp	_SuperposeCoordinates;
    core::HashTableEq_sp 	_Data;
	uint			_Members;

public:

	double	rmsCompareToEntry(Structure_Old_ListEntry_sp otherEntry);

	Structure_Old_List_sp	getStructureList() {_OF(); ASSERTNOTNULL(this->_WeakStructureList);return this->_WeakStructureList;};
	void	setStructureList(Structure_Old_List_sp s) { this->_WeakStructureList = s; };
	void	setMembers(uint i) { this->_Members = i; };
	uint	getMembers() { return this->_Members; };




	void	translateAllCoordinates(const Vector3& offset);

    core::HashTableEq_sp getData() { return this->_Data;};

        void setAllCoordinates(geom::CoordinateArray_sp ac);
        geom::CoordinateArray_sp getAllCoordinates() { return this->_AllCoordinates; }
        void setSuperposeCoordinates(geom::CoordinateArray_sp ac);
        geom::CoordinateArray_sp getSuperposeCoordinates() { return this->_SuperposeCoordinates; }

	void		applyAllCoordinatesToMatter(Matter_sp agg) { IMPLEMENT_ME();};
	void		applySuperposeCoordinatesToMatter(Matter_sp agg) { IMPLEMENT_ME();};

	void	writeCoordinatesToMatter(Matter_sp agg);
	void	extractCoordinatesFromMatter(Matter_sp agg);


	DEFAULT_CTOR_DTOR(Structure_Old_ListEntry_O);
};






class Structure_Old_List_O;
    typedef	gctools::smart_ptr<Structure_Old_List_O>	Structure_Old_List_sp;
SMART(StructureList);
class Structure_Old_List_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,Structure_Old_List_O,"StructureList",core::CxxObject_O);

public:
	void initialize();
public:
//	void	archiveBase(core::ArchiveP node);
public:
//	void	archive(core::ArchiveP node);
protected:
	Matter_sp				_Matter;
    gctools::Vec0<Structure_Old_ListEntry_sp>		_Entries;
    gctools::SmallOrderedSet<Atom_sp>			_AllAtoms;
    gctools::SmallOrderedSet<Atom_sp>			_SuperposeAtoms;
	double					_RmsCutOff;
    core::HashTableEq_sp			_Data;
protected:

public:
//    static Structure_Old_List_sp open_StructureList(const string& fn);

public:

    core::HashTableEq_sp	getData() { return this->_Data;};
    typedef	gctools::Vec0<Structure_Old_ListEntry_sp>::iterator	entryIterator;

	void	saveAs(const string& fn);


	void	setRmsCutOff(double co)	{this->_RmsCutOff = co;};

		/* Set the Matter object that this StructureList
		 * compares conformations of
		 */
	void	setMatter(Matter_sp agg);
	Matter_sp	getMatter() { return this->_Matter;};

		//! Clear the superposition atoms
	void	clearSuperposeAtoms() { this->_SuperposeAtoms.clear();};
		//! Add a superposition atom
    void	superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms);
		//! Set up to superpose all heavy atoms
	void	superposeAllHeavyAtoms();
	void	addSuperposeAtom(Atom_sp a);

    gctools::SmallOrderedSet<Atom_sp>::iterator begin_AllAtoms()
		{ return this->_AllAtoms.begin();};
    gctools::SmallOrderedSet<Atom_sp>::iterator end_AllAtoms()
		{ return this->_AllAtoms.end();};

    gctools::SmallOrderedSet<Atom_sp>::iterator begin_SuperposeAtoms()
		{ return this->_SuperposeAtoms.begin();};
    gctools::SmallOrderedSet<Atom_sp>::iterator end_SuperposeAtoms()
		{ return this->_SuperposeAtoms.end();};

    gctools::SmallOrderedSet<Atom_sp>	getSuperposeAtoms() { return this->_SuperposeAtoms; };

	entryIterator	begin_Entries()
		{ return this->_Entries.begin();};
	entryIterator	end_Entries()
		{ return this->_Entries.end();};

	uint	numberOfSuperposeAtoms()	{return this->_SuperposeAtoms.size();};
	virtual void	clearEntries();
	Structure_Old_ListEntry_sp getEntry(uint i);

	void	translateAllCoordinates(const Vector3& offset);

		/* Return true if the conformation is a new one and
		 * should be added to the StructureList
		 */

	Structure_Old_ListEntry_sp createStructureListEntryIfConformationIsNew(Matter_sp matter);

	/*! Add the entry to the StructureList, return status
	 * Subclass this to add the entry in a sorted order
	 * or to not add entries if they don't meet some criteria
	 */
    virtual uint addEntry(Structure_Old_ListEntry_sp entry);


	/*! Try to add each entry in structureList to this one.
	 * Since this will have a different SuperposeAtoms list we should
	 * get a subset of the structures.
	 */
	void	addAllEntries(Structure_Old_List_sp structureList);

	uint	numberOfEntries() { return this->_Entries.size();};


    virtual bool canRender() { return true; };
#ifdef RENDER
    virtual geom::Render_sp rendered(core::List_sp options);
#endif







    DEFAULT_CTOR_DTOR(Structure_Old_List_O);
};


};


TRANSLATE(chem::Structure_Old_ListEntry_O);
TRANSLATE(chem::Structure_Old_List_O);
#endif
