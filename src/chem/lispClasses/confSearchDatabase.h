       
       


/*
 */

#ifndef ConformationalSearchDatabase_H
#define	ConformationalSearchDatabase_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "addon/vector3.h"
#include "superposableConformationCollection.h"
//#include "bond.h"
//#include "atom.h"
//#include "residue.h"
//#include "nVector.h"
//#include "ffBaseDb.h"
//#include "ffTypesDb.h"
//#include "ffStretchDb.h"
//#include "ffAngleDb.h"
//#include "addon/coordinateArray.h"

#include "addon/quickDom.fwd.h"// confSearchDatabase.h wants QDomNode needs quickDom.fwd.h

namespace chem
{





typedef	enum	{ energyTooHigh=0,success=1} ConformationalSearchDatabaseFailReasonEnum;


SMART(QDomNode);
SMART(Atom);


SMART(ConformationalSearchDatabase);
class ConformationalSearchDatabase_O : public SuperposableConformationCollection_O
{
    LISP_BASE1(SuperposableConformationCollection_O);
    LISP_CLASS(chem,ChemPkg,ConformationalSearchDatabase_O,"ConformationalSearchDatabase");

public:
	void initialize();
public:
//	void	archive(core::ArchiveP node);
private:
	bool					_DisplayProgress;
        double					_EnergyCutOff;

private:
	void	_fillAtoms(addon::QDomNode_sp node, Aggregate_sp agg, VectorAtom& atoms);
	void	_fillEntries(addon::QDomNode_sp node);
	void	_writeEntryCoordinates(Aggregate_sp agg, int i);
public:

	void	enableDisplayProgress() { this->_DisplayProgress = true; };
	void	disableDisplayProgress() { this->_DisplayProgress = false; };
	void	setEnergyCutOff(double co) { this->_EnergyCutOff = co;};

	void	clearEntries();

	uint addEntry(ConformationCollectionEntry_sp entry);

	void		finishSearch();






	DEFAULT_CTOR_DTOR(ConformationalSearchDatabase_O);
};


};


TRANSLATE(chem::ConformationalSearchDatabase_O);
#endif
