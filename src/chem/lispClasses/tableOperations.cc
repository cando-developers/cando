/*
    File: tableOperations.cc
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
#define DEBUG_LEVEL_NONE

#include "core/common.h"
#include "core/symbol.h"
#include "core/evaluator.h"
#include "core/sort.h"
#include "loop.h"
#include "symbolTable.h"
#include "core/environment.h"
#include "core/activationFrame.h"
#include "tableOperations.h"
#include "superposeEngine.h"
#include "core/wrappers.h"

namespace chem
{

DOCGROUP(cando)
CL_DEFUN bool chem__verifyThatFieldOfAllEntriesAreEqual(Table_sp table, core::Symbol_sp field )
    {
	if ( table->numberOfEntries() <= 1 ) return true;
	TableEntry_sp firstEntry = *(table->beginEntry());
	core::T_sp firstObject = firstEntry->read(field);
	uint idx=1;
	for (Table_O::entryIterator it=table->beginEntry()+1; it!=table->endEntry(); it++,idx++ ) 
	{
	    core::T_sp otherObject = (*it)->read(field);
	    if ( !firstObject->equal(otherObject) ) 
	    {
              LOG(BF("The firstEntry field(%s)/value(%s) is not equal to that of entry(%d)/value(%s)") % field->fullName() % _rep_(firstObject) % idx % _rep_(otherObject) );
		return false;
	    }
	}
	LOG(BF("They all match"));
	return true;
    }



DOCGROUP(cando)
CL_DEFUN bool chem__superposeAllHeavyAtoms(Table_sp table, core::Symbol_sp matterField, core::Symbol_sp writeRmsd )
    {
	if ( !table->hasField(writeRmsd) )
	{
	    table->appendField(writeRmsd,core::DoubleFloat_O::___staticClass);
	}
	if ( !(table->fieldClass(writeRmsd) == core::DoubleFloat_O::___staticClass ))
	{
	    SIMPLE_ERROR(BF("The tables field for writing the RMSD is not Real!"));
	}
	if ( !verifyThatFieldOfAllEntriesAreEqual(table,matterField) )
	{
	    SIMPLE_ERROR(BF("If you want to superpose fields of a table they must all be equal and these are not"));
	}
	TableEntry_sp first = *(table->beginEntry());
    	//
	// The RMSD of the first entry of the table is always 0.0
	//
	core::DoubleFloat_sp val0 = core::DoubleFloat_O::create(0.0);
	first->write(writeRmsd,val0);
	if ( table->numberOfEntries() <= 1 ) return false;

	// 
	// Create the Superposer
	SuperposeSelectedAtoms_sp sm = SuperposeSelectedAtoms_O::create();
	sm->setMatter(first->read(matterField).as<Matter_O>());
	Loop lAtoms(sm->getMatter(),ATOMS);
	SYMBOL_EXPORT_SC_(ChemPkg,superpose);
	core::Symbol_sp superposeSymbol = _sym_superpose;
	while ( lAtoms.advance() )
	{
	    Atom_sp a = lAtoms.getAtom();
	    a->clearProperty(superposeSymbol);
	    if ( a->isHeavyAtom() )
	    {
		a->setPropertyTrue(superposeSymbol);
	    }
	}
	sm->updateSuperposeAtoms();
	sm->copyMatterCoordinatesIntoFixedCoordinates(first->read(matterField).as<Matter_O>());
	Table_O::entryIterator ei;
	for ( ei=table->beginEntry()+1; ei!=table->endEntry(); ei++ )
	{
	    sm->copyMatterCoordinatesIntoMoveableCoordinates((*ei)->read(matterField).as<Matter_O>());
	    sm->doSuperpose();
	    (*ei)->write(writeRmsd,core::DoubleFloat_O::create(sm->rootMeanSquareDifference()));
	    (*ei)->read(matterField).as<Matter_O>()->applyTransformToAtoms(sm->superpose());
	}
	return true;
    }






class	OrderUsingCallback
{
private:
	core::Function_sp	_compareForm;
public:
    OrderUsingCallback(core::Function_sp compareForm )
    {
	this->_compareForm = compareForm;
    }
    bool operator()(TableEntry_sp x, TableEntry_sp y )
    {
	core::T_sp result = core::eval::funcall(this->_compareForm,x,y);
	if ( result.isTrue() ) return true;
	return false;
    }
};


/*!
 * The comparison function must take two TableEntries and return true is the first entry
 * comes before the second one.
 */
    void sortTableEntriesUsingComparisonFunction(Table_sp table, core::Function_sp compare )
    {
	OrderUsingCallback orderer(compare);
	core::sort::quickSort(table->beginEntry(), table->endEntry(),orderer );
    }


/*!
 * The comparison function must take two TableEntries and return true is the first entry
 * comes before the second one.
 */
DOCGROUP(cando)
CL_DEFUN void chem__forEachEntryEvaluate(Table_sp table, core::Function_sp eval )
    {
    for ( Table_O::entryIterator ei=table->beginEntry(); ei!=table->endEntry(); ei++ )
    {
	core::Cons_sp form = core::Cons_O::createList(eval,*ei);
	core::eval::evaluate(form,_Nil<core::Environment_O>());
    }
}




void initialize_tableOperations(core::LispPtr lisp)
{
//    core::af_def( ChemPkg, "verifyThatFieldOfAllEntriesAreEqual", &verifyThatFieldOfAllEntriesAreEqual);
//    core::af_def( ChemPkg, "sortTableEntriesUsingComparisonFunction", &sortTableEntriesUsingComparisonFunction);
//    core::af_def( ChemPkg, "superposeAllHeavyAtoms", &superposeAllHeavyAtoms);
//    core::af_def( ChemPkg, "forEachEntryEvaluate", &forEachEntryEvaluate);
}

};
