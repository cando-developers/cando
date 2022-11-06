/*
    File: conformationExplorer.h
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

#ifndef ConformationExplorer_H
#define	ConformationExplorer_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/geom/vector3.h>
#include <clasp/core/environment.fwd.h>
#include <cando/geom/coordinateArray.fwd.h>

namespace       chem {

extern	string	ConformationExplorerEntry_Xml;




SMART(ConformationExplorerEntryStage);
SMART(ConformationExplorerEntry);
SMART(ConformationExplorerMatch);
SMART(ConformationExplorer);
SMART(SuperposeEngine);


#define	ConformationExplorer_AddedEntry	0


SMART(ConformationExplorerEntryStage);
class ConformationExplorerEntryStage_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,ConformationExplorerEntryStage_O,"ConformationExplorerEntryStage",core::CxxObject_O);

public:
	void	initialize();
private:
  core::T_sp		_BasedOnStageName;
  core::T_sp		_StageName;
  bool				_Complete;
  ConformationExplorerEntry_sp	_WeakConformationExplorerEntry;
  geom::SimpleVectorCoordinate_sp	_FinalCoordinates;
  core::HashTableEq_sp			_Binder;
			//! Energy of the conformation in kcal as calculated by the package/model
  double				_EnergyKCal;
			//! ExternalInterfaceName is cando, gamess, gaussian, nwchem etc
  core::T_sp		_ExternalInterfaceName;
			//! Model is gaff, am1, pm3, b3lyp etc.
  core::T_sp		_Model;
public:
        void setFinalCoordinates(geom::SimpleVectorCoordinate_sp ac);
CL_LISPIFY_NAME("getFinalCoordinates");
CL_DEFMETHOD         geom::SimpleVectorCoordinate_sp getFinalCoordinates() { return this->_FinalCoordinates; }
public:

CL_LISPIFY_NAME("getEnergyKCal");
CL_DEFMETHOD 	double		getEnergyKCal() { return this->_EnergyKCal; };
CL_LISPIFY_NAME("setEnergyKCal");
CL_DEFMETHOD 	void		setEnergyKCal(double energy) { this->_EnergyKCal = energy; };

CL_LISPIFY_NAME("getExternalInterfaceName");
  CL_DEFMETHOD 	core::T_sp	getExternalInterfaceName() { return this->_ExternalInterfaceName;};
CL_LISPIFY_NAME("setExternalInterfaceName");
  CL_DEFMETHOD 	void		setExternalInterfaceName(core::T_sp package) { this->_ExternalInterfaceName = package;};

CL_LISPIFY_NAME("getModel");
  CL_DEFMETHOD 	core::T_sp	getModel() { return this->_Model;};
CL_LISPIFY_NAME("setModel");
  CL_DEFMETHOD 	void		setModel(core::T_sp package) { this->_Model = package;};


	bool	recognizesAtom(Atom_sp atom);
	void	setCoordinateForAtom(Atom_sp atom, Vector3 pos );

	ConformationExplorerEntry_sp	getConformationExplorerEntry();
	ConformationExplorer_sp	getConformationExplorer();

CL_LISPIFY_NAME("getBasedOnStageName");
 CL_DEFMETHOD 	core::T_sp	getBasedOnStageName() { return this->_BasedOnStageName; };
CL_LISPIFY_NAME("setBasedOnStageName");
 CL_DEFMETHOD 	void	setBasedOnStageName(core::T_sp key ) { this->_BasedOnStageName = key; };

CL_LISPIFY_NAME("getStageName");
 CL_DEFMETHOD 	core::T_sp	getStageName() { return this->_StageName; };
CL_LISPIFY_NAME("setStageName");
 CL_DEFMETHOD 	void	setStageName(core::T_sp key ) { this->_StageName = key; };

CL_LISPIFY_NAME("isComplete");
CL_DEFMETHOD 	bool	isComplete() { return this->_Complete; };
CL_LISPIFY_NAME("setComplete");
CL_DEFMETHOD 	void	setComplete(bool b) { this->_Complete = b;};

	void	setConformationExplorerEntry(ConformationExplorerEntry_sp s);
CL_LISPIFY_NAME("getBinder");
CL_DEFMETHOD 	core::HashTableEq_sp getBinder() { return this->_Binder; };

	void	translateFinalCoordinates(const Vector3& offset);

	void	writeCoordinatesToMatter(Matter_sp matter);
	void	extractCoordinatesFromMatter(Matter_sp matter);

	void	setSuperposableCoordinatesAsFixedWithinSuperposeEngine(SuperposeEngine_sp sup);
	void	setSuperposableCoordinatesAsMoveableWithinSuperposeEngine(SuperposeEngine_sp sup);

	void	setFinalCoordinatesAsFixedWithinSuperposeEngine(SuperposeEngine_sp sup);
	void	setFinalCoordinatesAsMoveableWithinSuperposeEngine(SuperposeEngine_sp sup);

#if 0
	bool	canRender() { return true;};
#ifdef RENDER
	geom::Render_sp	rendered(core::List_sp options);
#endif
#endif
	virtual	string	coordinatesAsString()	{_OF();SUBCLASS_MUST_IMPLEMENT();};
	virtual	bool	parseStructureFromFileName(const string& fileName) {_OF();SUBCLASS_MUST_IMPLEMENT();};

		/*! Align the coordinates to the given entryStage
		 * using the Superpose atoms of the ConformationExplorer
		 */
	void alignCoordinatesTo(ConformationExplorerEntryStage_sp entryStage );




	DEFAULT_CTOR_DTOR(ConformationExplorerEntryStage_O);
};








SMART(ConformationExplorerEntry);
class ConformationExplorerEntry_O : public core::CxxObject_O
{
    friend class ConformationExplorer_O;
    LISP_CLASS(chem,ChemPkg,ConformationExplorerEntry_O,"ConformationExplorerEntry",core::CxxObject_O);
public:
	void	initialize();
private:
	ConformationExplorer_sp				_WeakConformationExplorer;
    gctools::Vec0<ConformationExplorerEntryStage_sp>		_Stages;
		/*!One stage in the list of stages is selected
		 * for extraction of MonomerCoordinates or for writing out
		 * coordinates for external packages like Gaussian and Gamess to further optimize
		 */
	ConformationExplorerEntryStage_sp		_SelectedStage;
	/*! Every time an entry is created it gets a new UniqueEntryIndex */
	int							_UniqueEntryIndex;
	core::HashTableEq_sp				_Binder;
public:
    typedef	gctools::Vec0<ConformationExplorerEntryStage_sp>::iterator	stageIterator;
public:
CL_LISPIFY_NAME("getConformationExplorer");
CL_DEFMETHOD 	ConformationExplorer_sp	getConformationExplorer() {_OF(); ASSERTNOTNULL(this->_WeakConformationExplorer);return this->_WeakConformationExplorer;};
	void	setConformationExplorer(ConformationExplorer_sp s);
	core::HashTableEq_sp getBinder() { return this->_Binder; };
public:

	stageIterator	begin_Stages() { return this->_Stages.begin(); };
	stageIterator	end_Stages() { return this->_Stages.end(); };


		/*! Return true if there is a selected stage */
CL_LISPIFY_NAME("hasSelectedStage");
CL_DEFMETHOD 	bool	hasSelectedStage() {_OF(); ASSERTNOTNULL(this->_SelectedStage); return this->_SelectedStage.notnilp(); };
		/*! Return the selected stage, guaranteed to be one of this entries stages*/
	ConformationExplorerEntryStage_sp	getSelectedStage();
		/*! Set the selected stage, throws error if its not one of this entries stages */
	void	setSelectedStage(ConformationExplorerEntryStage_sp stage);

        core::HashTable_sp	getEntryStageNames();

CL_LISPIFY_NAME("getUniqueEntryIndex");
CL_DEFMETHOD 	int	getUniqueEntryIndex() const { return this->_UniqueEntryIndex;};

 core::T_sp			lastEntryStageName();
 bool			hasEntryStageWithName(core::T_sp key );
 ConformationExplorerEntryStage_sp	getEntryStage(core::T_sp key);
 
		/*! Return true if the last stage has the stageName and
		 * is complete
		 */
 bool	alreadyHasLastCompleteStage(core::T_sp stageName);

		/*! If the last EntryStage is incomplete then return it
		 * Otherwise create a new one and return it
		 * Give it the stageName and set it to incomplete
		 */
 ConformationExplorerEntryStage_sp	getOrCreateLastIncompleteEntryStage(core::T_sp stageName);

	ConformationExplorerEntryStage_sp	getLastCompleteEntryStage();
	ConformationExplorerEntryStage_sp	getLastEntryStage();

#if 0
	bool	canRender() { return true;};
#ifdef RENDER
	geom::Render_sp	rendered(core::List_sp options);
#endif
#endif
	virtual ConformationExplorerEntryStage_sp createEntryStage(core::T_sp name);



	DEFAULT_CTOR_DTOR(ConformationExplorerEntry_O);
};









SMART(ConformationExplorer);
class ConformationExplorer_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,ConformationExplorer_O,"ConformationExplorer",core::CxxObject_O);
public:
	void initialize();
public:
friend class ConformationExplorerEntry_O;
friend class ConformationExplorerEntryStage_O;
protected:
	Matter_sp				_Matter;
    gctools::Vec0<ConformationExplorerEntry_sp>	_Entries;
    gctools::SmallOrderedSet<Atom_sp>			_AllAtoms;
	core::ComplexVector_byte32_t_sp		_SuperposeAtomIndices;
	core::HashTableEq_sp			_Binder;
protected:
	Atom_sp	_getAtomAtIndex(unsigned i);
public:
    static ConformationExplorer_sp open_ConformationExplorer(const string& fn);

public:
    typedef	gctools::Vec0<ConformationExplorerEntry_sp>::iterator	entryIterator;
    typedef	gctools::SmallOrderedSet<Atom_sp>::iterator			atomIterator;

public:
	geom::SimpleVectorCoordinate_sp	_SimpleVectorCoordinate(Matter_sp matter);
	core::ComplexVector_byte32_t_sp		_getSuperposeAtomIndices() { return this->_SuperposeAtomIndices;};
public:

		//! Clear the superposition atoms
	void	clearSuperposeAtoms();
		//! Add a superposition atom
    void	superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms);
		//! Set up to superpose all heavy atoms
	void	superposeAllHeavyAtoms();
	void	addSuperposeAtom(Atom_sp a);

    gctools::SmallOrderedSet<Atom_sp>	getSuperposeAtoms();
	uint	numberOfSuperposeAtoms();

public:
CL_LISPIFY_NAME("getBinder");
CL_DEFMETHOD 	core::HashTableEq_sp	getBinder() { return this->_Binder;};

	void	saveAs(const string& fn);

		/*! Set the Matter object that this ConformationExplorer
		 * compares conformations of
		 */
	virtual void	setMatter(Matter_sp agg);
CL_LISPIFY_NAME("getMatter");
CL_DEFMETHOD 	Matter_sp	getMatter() { return this->_Matter;};

	atomIterator begin_AllAtoms()
		{ return this->_AllAtoms.begin();};
	atomIterator end_AllAtoms()
		{ return this->_AllAtoms.end();};

	entryIterator	begin_Entries()
		{ return this->_Entries.begin();};
	entryIterator	end_Entries()
		{ return this->_Entries.end();};

		//! Return the index of the entry
	unsigned getEntryIndex(ConformationExplorerEntry_sp entry);

	core::List_sp entriesAsList();

	virtual void	clearEntries();
	ConformationExplorerEntry_sp getEntry(unsigned i);
	ConformationExplorerEntry_sp firstEntry();

	ConformationExplorerEntry_sp getClosestMatchingEntry(ConformationExplorerMatch_sp match);
	ConformationExplorerEntryStage_sp getClosestMatchingEntryStage(ConformationExplorerMatch_sp match);

	void	translateFinalCoordinates(const Vector3& offset);

		/* Return true if the conformation is a new one and
		 * should be added to the ConformationExplorer
		 */

	/*! Add the entry to the ConformationExplorer, return status
	 * Subclass this to add the entry in a sorted order
	 * or to not add entries if they don't meet some criteria
	 */
    	void appendEntry(ConformationExplorerEntry_sp entry);

	bool	hasStageNameInAllEntries(core::T_sp stageName);
	bool	findMostSimilarConformationEntryStageWithStageName(
			Matter_sp 			matter,
			core::T_sp 			stageName,
			double&				rmsDifference,
			ConformationExplorerEntryStage_sp &closestEntryStage,
			uint&				bestEntryIndex );

	void	findClosestMatchingConformation(ConformationExplorerMatch_sp match, Matter_sp matter, core::T_sp stageName, double rms);


	/*! Sort the ConformationExplorer based on a value in one of the entry stages
	 */
	void	sort(core::T_sp stageName, core::Symbol_sp keyValue);


//	void	addConformationExplorer(ConformationExplorer_sp structureList);

CL_LISPIFY_NAME("numberOfEntries");
CL_DEFMETHOD 	uint numberOfEntries() { return this->_Entries.size();};
	uint numberOfAllAtoms() { return this->_AllAtoms.size();};

	bool	recognizesAtomName(core::T_sp name);
	uint indexOfAtom(Atom_sp atom);

	void	alignAllConformationsToTheFirstForStage(core::T_sp stageName );

#if 0
	bool canRender() { return true; }
#ifdef RENDER
	geom::Render_sp rendered(core::List_sp options);
#endif
#endif

		/*! Create a new entry but dont add it yet
		 * Subclasses of ConformationExplorer can return
		 * subclasses of ConformationExplorerEntry with different
		 * behavior
		 */
virtual	ConformationExplorerEntry_sp	createEntry();

		/*! Sort the entries by the energy of the given stage
		 */
 void sortByEnergyAscendingForStage(core::T_sp stageName);






	DEFAULT_CTOR_DTOR(ConformationExplorer_O);
};



/*!
 * When you search a ConformationExplorer this stores the index of the conformation
 * that it most closely matches and whether it matches within the rmsd specified
 */
SMART(ConformationExplorerMatch );
class ConformationExplorerMatch_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,ConformationExplorerMatch_O,"ConformationExplorerMatch",core::CxxObject_O);
public: // virtual functions inherited from Object
	void	initialize();
//	string	__repr__() const;

private: // instance variables
	bool			_Matches;
	uint			_ClosestMatchIndex;
	double			_ClosestMatchRms;
  core::T_sp	_MatchStageName;
public:	// Creation class functions

public:

	void setMatches(bool b) { this->_Matches = b; };
CL_LISPIFY_NAME("getMatches");
CL_DEFMETHOD 	bool getMatches() { return this->_Matches; };
 void setMatchStageName(core::T_sp s) { this->_MatchStageName=s;};
 core::T_sp getMatchStageName() { return this->_MatchStageName;};
CL_LISPIFY_NAME("getClosestMatchRms");
CL_DEFMETHOD 	double getClosestMatchRms() { return this->_ClosestMatchRms;};
	void setClosestMatchRms(double i) { this->_ClosestMatchRms= i;};
CL_LISPIFY_NAME("getClosestMatchIndex");
CL_DEFMETHOD 	uint getClosestMatchIndex() { return this->_ClosestMatchIndex;};
	void setClosestMatchIndex(uint i) { this->_ClosestMatchIndex = i;};
	ConformationExplorerMatch_O( const ConformationExplorerMatch_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ConformationExplorerMatch_O);
};








};


TRANSLATE(chem::ConformationExplorerEntryStage_O);
TRANSLATE(chem::ConformationExplorerEntry_O);
TRANSLATE(chem::ConformationExplorer_O);
TRANSLATE(chem::ConformationExplorerMatch_O);
#endif
