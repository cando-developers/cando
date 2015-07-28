       
       
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
#include <clasp/core/intArray.h>

#include <clasp/core/environment.fwd.h>
#include <cando/geom/coordinateArray.fwd.h>
#include <clasp/core/intArray.fwd.h>

namespace       chem {

extern	string	ConformationExplorerEntry_Xml;




SMART(ConformationExplorerEntryStage);
SMART(ConformationExplorerEntry);
SMART(ConformationExplorerMatch);
SMART(ConformationExplorer);
SMART(SuperposeEngine);
SMART(CoordinateArray);
SMART(IntArray);



#define	ConformationExplorer_AddedEntry	0


SMART(ConformationExplorerEntryStage);
class ConformationExplorerEntryStage_O : public core::CxxObject_O
{
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(chem,ChemPkg,ConformationExplorerEntryStage_O,"ConformationExplorerEntryStage");

public:
//    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
	string				_BasedOnStageName;
	string				_StageName;
	bool				_Complete;
	ConformationExplorerEntry_wp	_WeakConformationExplorerEntry;
        geom::CoordinateArray_sp	_FinalCoordinates;
	core::HashTableEq_sp			_Binder;
			//! Energy of the conformation in kCal as calculated by the package/model
	double				_EnergyKCal;
			//! ExternalInterfaceName is cando, gamess, gaussian, nwchem etc
	string				_ExternalInterfaceName;
			//! Model is gaff, am1, pm3, b3lyp etc.
	string				_Model;
public:
        void setFinalCoordinates(geom::CoordinateArray_sp ac);
        geom::CoordinateArray_sp getFinalCoordinates() { return this->_FinalCoordinates; }
public:

	double		getEnergyKCal() { return this->_EnergyKCal; };
	void		setEnergyKCal(double energy) { this->_EnergyKCal = energy; };

	string		getExternalInterfaceName() { return this->_ExternalInterfaceName;};
	void		setExternalInterfaceName(const string& package) { this->_ExternalInterfaceName = package;};

	string		getModel() { return this->_Model;};
	void		setModel(const string& package) { this->_Model = package;};


	bool	recognizesAtom(Atom_sp atom);
	void	setCoordinateForAtom(Atom_sp atom, Vector3 pos );

	ConformationExplorerEntry_sp	getConformationExplorerEntry();
	ConformationExplorer_sp	getConformationExplorer();

	string	getBasedOnStageName() { return this->_BasedOnStageName; };
	void	setBasedOnStageName(const string& key ) { this->_BasedOnStageName = key; };

	string	getStageName() { return this->_StageName; };
	void	setStageName(const string& key ) { this->_StageName = key; };

	bool	isComplete() { return this->_Complete; };
	void	setComplete(bool b) { this->_Complete = b;};

	void	setConformationExplorerEntry(ConformationExplorerEntry_sp s);
	core::HashTableEq_sp getBinder() { return this->_Binder; };

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
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(chem,ChemPkg,ConformationExplorerEntry_O,"ConformationExplorerEntry");
public:
//    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
	ConformationExplorer_wp				_WeakConformationExplorer;
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
	ConformationExplorer_sp	getConformationExplorer() {_OF(); ASSERTNOTNULL(this->_WeakConformationExplorer);return this->_WeakConformationExplorer;};
	void	setConformationExplorer(ConformationExplorer_sp s);
	core::HashTableEq_sp getBinder() { return this->_Binder; };
public:

	stageIterator	begin_Stages() { return this->_Stages.begin(); };
	stageIterator	end_Stages() { return this->_Stages.end(); };


		/*! Return true if there is a selected stage */
	bool	hasSelectedStage() {_OF(); ASSERTNOTNULL(this->_SelectedStage); return this->_SelectedStage.notnilp(); };
		/*! Return the selected stage, guaranteed to be one of this entries stages*/
	ConformationExplorerEntryStage_sp	getSelectedStage();
		/*! Set the selected stage, throws error if its not one of this entries stages */
	void	setSelectedStage(ConformationExplorerEntryStage_sp stage);

	adapt::StringSet_sp	getEntryStageNames();

	int	getUniqueEntryIndex() const { return this->_UniqueEntryIndex;};

	string			lastEntryStageName();
	bool			hasEntryStageWithName(const string& key );
	ConformationExplorerEntryStage_sp	getEntryStage(const string& key);

		/*! Return true if the last stage has the stageName and
		 * is complete
		 */
	bool	alreadyHasLastCompleteStage(const string& stageName);

		/*! If the last EntryStage is incomplete then return it
		 * Otherwise create a new one and return it
		 * Give it the stageName and set it to incomplete
		 */
	ConformationExplorerEntryStage_sp	getOrCreateLastIncompleteEntryStage(const string& stageName);

	ConformationExplorerEntryStage_sp	getLastCompleteEntryStage();
	ConformationExplorerEntryStage_sp	getLastEntryStage();

#if 0
	bool	canRender() { return true;};
#ifdef RENDER
	geom::Render_sp	rendered(core::List_sp options);
#endif
#endif
	virtual ConformationExplorerEntryStage_sp createEntryStage(string const& name);



	DEFAULT_CTOR_DTOR(ConformationExplorerEntry_O);
};









SMART(ConformationExplorer);
class ConformationExplorer_O : public core::CxxObject_O
{
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(chem,ChemPkg,ConformationExplorer_O,"ConformationExplorer");
public:
	void initialize();
public:
//	void	archiveBase(core::ArchiveP node);
friend class ConformationExplorerEntry_O;
friend class ConformationExplorerEntryStage_O;
protected:
	Matter_sp				_Matter;
    gctools::Vec0<ConformationExplorerEntry_sp>	_Entries;
    gctools::SmallOrderedSet<Atom_sp>			_AllAtoms;
	core::IntArray_sp				_SuperposeAtomIndices;
	core::HashTableEq_sp			_Binder;
protected:
	Atom_sp	_getAtomAtIndex(unsigned i);
public:
    static ConformationExplorer_sp open_ConformationExplorer(const string& fn);

public:
    typedef	gctools::Vec0<ConformationExplorerEntry_sp>::iterator	entryIterator;
    typedef	gctools::SmallOrderedSet<Atom_sp>::iterator			atomIterator;
	typedef	core::IntArray_O::iterator				superposeAtomIndexIterator;

public:
	geom::CoordinateArray_sp	_extractCoordinateArray(Matter_sp matter);
	core::IntArray_sp		_getSuperposeAtomIndices() { return this->_SuperposeAtomIndices;};
public:

	superposeAtomIndexIterator	superposeAtomIndices_begin() { return this->_SuperposeAtomIndices->begin(); };
	superposeAtomIndexIterator	superposeAtomIndices_end() { return this->_SuperposeAtomIndices->end(); };

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
	core::HashTableEq_sp	getBinder() { return this->_Binder;};

	void	saveAs(const string& fn);

		/*! Set the Matter object that this ConformationExplorer
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
	unsigned getEntryIndex(ConformationExplorerEntry_sp entry);

	core::List_sp entriesAsCons();

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

	bool	hasStageNameInAllEntries(const string& stageName);
	bool	findMostSimilarConformationEntryStageWithStageName(
			Matter_sp 			matter,
			const string& 			stageName,
			double&				rmsDifference,
			ConformationExplorerEntryStage_sp &closestEntryStage,
			uint&				bestEntryIndex );

	void	findClosestMatchingConformation(ConformationExplorerMatch_sp match, Matter_sp matter, const string& stageName, double rms);



#ifdef	USEBOOSTPYTHON
	boost::python::tuple boost_findMostSimilarConformationEntryStageWithStageName(Matter_sp matter, const string& stageName );
#endif

	/*! Sort the ConformationExplorer based on a value in one of the entry stages
	 */
	void	sort(const string& stageName, core::Symbol_sp keyValue);


//	void	addConformationExplorer(ConformationExplorer_sp structureList);

	uint numberOfEntries() { return this->_Entries.size();};
	uint numberOfAllAtoms() { return this->_AllAtoms.size();};

	bool	recognizesAtomName(const string& name);
	uint indexOfAtom(Atom_sp atom);

	void	alignAllConformationsToTheFirstForStage(const string& stageName );

#if 0
	bool canRender() { return true; }
#ifdef RENDER
	geom::Render_sp rendered(core::List_sp options);
#endif
#endif
	string __str__() ;

		/*! Create a new entry but dont add it yet
		 * Subclasses of ConformationExplorer can return
		 * subclasses of ConformationExplorerEntry with different
		 * behavior
		 */
virtual	ConformationExplorerEntry_sp	createEntry();

		/*! Sort the entries by the energy of the given stage
		 */
	void sortByEnergyAscendingForStage(const string& stageName);






	DEFAULT_CTOR_DTOR(ConformationExplorer_O);
};



/*!
 * When you search a ConformationExplorer this stores the index of the conformation
 * that it most closely matches and whether it matches within the rmsd specified
 */
SMART(ConformationExplorerMatch );
class ConformationExplorerMatch_O : public core::CxxObject_O
{
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(chem,ChemPkg,ConformationExplorerMatch_O,"ConformationExplorerMatch");
public: // virtual functions inherited from Object
	void	initialize();
//	string	__repr__() const;

private: // instance variables
	bool			_Matches;
	uint			_ClosestMatchIndex;
	double			_ClosestMatchRms;
	string			_MatchStageName;
public:	// Creation class functions

public:

	void setMatches(bool b) { this->_Matches = b; };
	bool getMatches() { return this->_Matches; };
	void setMatchStageName(const string& s) { this->_MatchStageName=s;};
	string getMatchStageName() { return this->_MatchStageName;};
	double getClosestMatchRms() { return this->_ClosestMatchRms;};
	void setClosestMatchRms(double i) { this->_ClosestMatchRms= i;};
	uint getClosestMatchIndex() { return this->_ClosestMatchIndex;};
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
