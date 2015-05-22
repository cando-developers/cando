       
       



#ifndef	MonomerCoordinates_H //[
#define MonomerCoordinates_H


#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "monomerContext.h"
#include "conformationExplorer.h"



#include "chemPackage.h"


namespace chem {

    FORWARD(Rotamer);
SMART(Topology);

SMART(MonomerCoordinates);
class MonomerCoordinates_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,MonomerCoordinates_O,"MonomerCoordinates");
#if INIT_TO_FACTORIES
 public:
    static MonomerCoordinates_sp make(MonomerContext_sp context);
#else
    DECLARE_INIT();
#endif
public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);

public:
    typedef	gctools::Vec0<Rotamer_sp>	Rotamers;
private:
	string				_Comment;
	MonomerContext_sp		_Context;
	string				_TrainerName;
	Rotamers			_Rotamers;
public:
	void _defineFromConformationExplorerOrDebug( Topology_sp topology,
				ConformationExplorer_sp structureList,
				Monomer_sp  focusMonomer,
				bool	   debug,
				uint maximumNumberOfConformations );
public:

	void setComment(const string& context) { this->_Comment = context; };
	string getComment() { return this->_Comment; };

	void setContext(MonomerContext_sp context) { this->_Context = context; };
	MonomerContext_sp getContext() { return this->_Context; };

	void defineFromConformationExplorer( Topology_sp topology,
				ConformationExplorer_sp structureList,
				Monomer_sp  focusMonomer,
				uint maximumNumberOfConformations );

	/*! Test if the extraction process will work with the current Aggregate
	 * if it fails then there is a problem in the MetaDatabase that needs to
	 * be fixed before we go forward
	 * Will return a Cons with two KeyedObjects (sawProblem: (bool))
	 * and (problems: (StringStream))
	 */

	 core::Cons_sp testExtraction( Topology_sp topology, 
			Aggregate_sp aggregate, 
			Monomer_sp focusMonomer );
	
	DEFAULT_CTOR_DTOR(MonomerCoordinates_O);
};




};
TRANSLATE(chem::MonomerCoordinates_O);
#endif //]


