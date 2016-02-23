/*
    File: confSearchDatabase.h
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

#ifndef ConformationalSearchDatabase_H
#define	ConformationalSearchDatabase_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "candoBase/vector3.h"
#include "superposableConformationCollection.h"
//#include "bond.h"
//#include "atom.h"
//#include "residue.h"
//#include "nVector.h"
//#include "ffBaseDb.h"
//#include "ffTypesDb.h"
//#include "ffStretchDb.h"
//#include "ffAngleDb.h"
//#include "candoBase/coordinateArray.h"

#include "candoBase/quickDom.fwd.h"// confSearchDatabase.h wants QDomNode needs quickDom.fwd.h

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
	void	_fillAtoms(candoBase::QDomNode_sp node, Aggregate_sp agg, VectorAtom& atoms);
	void	_fillEntries(candoBase::QDomNode_sp node);
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
