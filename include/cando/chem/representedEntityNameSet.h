#ifndef	RepresentedEntityNameSet_H
#define RepresentedEntityNameSet_H





#include <cando/chem/entityNameSet.h>


#include <cando/chem/chemPackage.h>


namespace chem {

SMART(MapOfMonomerNamesToAtomIndexers);

SMART(RepresentedEntityNameSet);
class RepresentedEntityNameSet_O : public EntityNameSet_O
{
    LISP_CLASS(chem,ChemPkg,RepresentedEntityNameSet_O,"RepresentedEntityNameSet",EntityNameSet_O);
#if INIT_TO_FACTORIES
 public:
    static RepresentedEntityNameSet_sp make(core::Symbol_sp representativeName);
#else
public:
    DECLARE_INIT();
#endif
public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);

private:
	core::Symbol_sp	_Representative;
public:
#if 0
    static RepresentedEntityNameSet_sp	create(CandoDatabase_sp db)
    {
        GC_ALLOCATE(RepresentedEntityNameSet_O, group ); // RP_Create<RepresentedEntityNameSet_O>(e);
//	group->setCandoDatabase(db);
	return group;
    };
#endif
public:
    virtual RepresentativeList_sp expandedRepresentativeList() const;

public:
    string __repr__() const { return this->description(); };
	string description() const;
	bool	hasRepresentative();
	void	setRepresentative(core::Symbol_sp cc );
	core::Symbol_sp	getRepresentative();

	/*! Return a copy of this object */
	EntityNameSetBase_sp copy() const;


		/*! Return a RepresentedEntityNameSet containing itself as its
		 * contents - this is used when creating trainers
		 * so that the contents of the RepresentedEntityNameSet are not locked
		 * in at the time the trainer is defined.
		 */
//	RepresentedEntityNameSet_sp asUnexpandedRepresentedEntityNameSet();

CL_NAME("testRepresentedEntityNameSet");
CL_DEFMETHOD 	string	testRepresentedEntityNameSet() { return "testRepresentedEntityNameSet"; };
//	virtual	string	getGeneralKey();

	RepresentedEntityNameSet_O( const RepresentedEntityNameSet_O& ss ); //!< Copy constructor



	DEFAULT_CTOR_DTOR(RepresentedEntityNameSet_O);
};






};


TRANSLATE(chem::RepresentedEntityNameSet_O);
#endif
