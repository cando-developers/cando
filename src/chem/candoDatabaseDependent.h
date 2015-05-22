       
       
#ifndef	CandoDatabaseDependent_H //[
#define CandoDatabaseDependent_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "core/symbolSet.fwd.h"



#include "chemPackage.h"


namespace chem {

SMART(CandoDatabase);
SMART(Constitution);

#if 0
SMART(CandoDatabaseDependent);
__ BEGIN_CLASS_DEFINITION(CandoDatabase_ODependent,CandoDatabaseDependent,O EnvironmentDependent) // {
public:
	void initialize();
public:
//	void	archiveBase(core::ArchiveP node);
//	bool	loadFinalize(core::ArchiveP node);
private:
//	CandoDatabase_wp	_WeakCandoDatabase;
public:


//    virtual void	setCandoDatabase(CandoDatabase_sp bdb);

		/*! Return true if this name can not be expanded.
		 * For instance monomer groups can be expanded to every monomer within them.
		 */
    virtual bool	isTerminalName() { return false; };
    		/*! Get the name of this object
		 */
    virtual string	getName()	{_OF();SUBCLASS_MUST_IMPLEMENT();};

    		/*! Return a list of names that is expanded recursively as much as possible
		 */
    virtual core::SymbolSet_sp	expandedNameSet() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
    		/*! Return a list of names that is expanded recursively but terminate with
		 * representative names.
		 */
    virtual RepresentativeList_sp expandedRepresentativeList() const {_OF(); SUBCLASS_MUST_IMPLEMENT(); };

		/*! If this object is a Stereoisomer then it has a Constitution otherwise it doesn't */
    virtual bool hasConstitution() { return false; };
    		/*! Return the constitution */
    virtual Constitution_sp constitution();


	CandoDatabase_sp getCandoDatabase();
	bool hasCandoDatabase();


	CandoDatabase_ODependent( const CandoDatabase_ODependent& ss ); //!< Copy constructor


__END_CLASS_DEFINITION(CandoDatabase_ODependent) //}
#endif

};
#endif //]
