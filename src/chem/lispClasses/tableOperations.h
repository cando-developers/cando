/*
    File: tableOperations.h
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

#ifndef	tableOperations_H
#define tableOperations_H



#include "chemPackage.h"
#include "core/foundation.h"
#include "core/lisp.fwd.h"// tableOperations.h wants Lisp needs lisp.fwd.h
#include "core/symbol.fwd.h"// tableOperations.h wants Symbol needs symbol.fwd.h
#include "core/executables.fwd.h"// tableOperations.h wants Executable needs intrinsics.fwd.h

namespace chem
{
SMART(Table);
SMART(TableEntry);
SMART(Lisp);
SMART(Symbol);
SMART(Executable);

	/*! For each entry in the table check to make sure that 
	 * the value in the given field is "equal()" to the value
	 * in that field of the first value.
	 * This allows you to verify that every entry has the same molecule
	 * in a field.
	 */
    bool verifyThatFieldOfAllEntriesAreEqual(Table_sp table, core::Symbol_sp field );

	/*! Superpose all molecules (using heavy atoms only) to the first entry
	 * and write the RMSD to the field (writeRmsd).
	 * If writeRmsd is nil then don't write the RMSD
	 * Throw an exception if all of the Matters in (superposeField) are not equal.
	 */
bool superposeAllHeavyAtoms(Table_sp table, core::Symbol_sp superposeField, core::Symbol_sp writeRmsd );

	/*! Evaluate the passed function for each entry
	 */
void forEachEntryEvaluate(Table_sp table, core::Function_sp eval );


	/*! Sort the table entries using a comparison function
	 * that accepts two TableEntry arguments and returns true if the
	 * first should come before the second
	 */
void sortTableEntriesUsingComparisonFunction(Table_sp table, core::Function_sp compare );

	/*! Sort the table entries using fields.
	 * Provide a list of fields and sort order in a list of two element lists.
	 * orderBy '( ( :energy :ascending ) (:rmsd :descending) )
	 * in ascending or descending order.
	 */
void sortTableEntriesUsingComparisonFunction(Table_sp table, core::Function_sp compare );


void initialize_tableOperations(core::LispPtr);

};
#endif
