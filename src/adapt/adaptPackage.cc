/*
    File: adaptPackage.cc
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
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/package.h>
#include <clasp/core/environment.h>
#include <cando/adapt/adaptPackage.h>
#include <cando/adapt/indexedObjectBag.h>
#include <cando/adapt/objectSet.h>
#include <cando/adapt/stringList.h>
#include <cando/adapt/symbolSet.h>
#include <clasp/core/object.h>
#include <clasp/asttooling/astVisitor.h>
//#include <clasp/cffi/cffi.h>
#include <clasp/clbind/class_registry.h>
#include <clasp/core/serialize.h>
#include <clasp/core/conditions.h>
#include <clasp/core/character.h>
#include <clasp/core/cons.h>
#include <clasp/core/cxxObject.h>
#include <clasp/core/fileSystem.h>
#include <clasp/core/environment.h>
#include <clasp/core/externalObject.h>
#include <clasp/core/hashTable.h>
#include <clasp/core/intArray.h>
#include <clasp/core/primitives.h>
#include <clasp/core/iterator.h>
#include <clasp/core/lambdaListHandler.h>
#include <clasp/core/loadTimeValues.h>
#include <clasp/core/multiStringBuffer.h>
#include <clasp/core/numbers.h>
#include <clasp/core/package.h>
#include <clasp/core/pathname.h>
#include <clasp/core/pointer.h>
#include <clasp/core/posixTime.h>
#include <clasp/core/random.h>
#include <clasp/core/readtable.h>
#include <clasp/core/reader.h>
#include <clasp/core/record.h>
#include <clasp/core/singleDispatchMethod.h>
#include <clasp/core/smallMap.h>
#include <clasp/core/smallMultimap.h>
#include <clasp/core/sourceFileInfo.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolToEnumConverter.h>
#include <clasp/core/symbol.h>
#include <clasp/core/lispList.h>
#include <clasp/core/weakHashTable.h>
#include <clasp/core/weakKeyMapping.h>
#include <clasp/core/weakPointer.h>
#include <clasp/core/wrappedPointer.h>
#include <clasp/llvmo/llvmoExpose.h>
#include <clasp/llvmo/debugInfoExpose.h>
#include <clasp/llvmo/debugLoc.h>
#include <clasp/llvmo/insertPoint.h>
#include <cando/adapt/iterateCons.h>
#include <cando/adapt/quickDom.h>
#include <cando/adapt/stringSet.h>
#include <cando/adapt/symbolList.h>
#include <cando/chem/ringFinder.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/alias.h>
#include <cando/chem/angle.h>
#include <cando/chem/atomIdMap.h>
#include <cando/chem/atomId.h>
#include <cando/chem/atomIndexer.h>
#include <cando/chem/atomReference.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/bond.h>
#include <cando/chem/chemdraw.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/chemInfo.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/command.h>
#include <cando/chem/complexRestraints.h>
#include <cando/chem/conformationCollection.h>
#include <cando/chem/conformationExplorer.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/coordSys.h>
#include <cando/chem/coupling.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/entity.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/frameRecognizer.h>
#include <cando/chem/improperTorsion.h>
#include <cando/chem/iterateMatter.h>
#include <cando/chem/iterateRestraints.h>
#include <cando/chem/macroModel.h>
#include <cando/chem/matter.h>
#include <cando/chem/minimizerLog.h>
#include <cando/chem/minimizer.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/monomerCoordinates.h>
#include <cando/chem/nVector.h>
#include <cando/chem/numericalFunction.h>
#include <cando/chem/candoScript.h>
#include <cando/chem/oligomer.h>
#include <cando/chem/monomer.h>
#include <cando/chem/typeAssignmentRules.h>
#include <cando/chem/pdbMonomer.h>
#include <cando/chem/pdb.h>
#include <cando/chem/plug.h>
#include <cando/chem/properTorsion.h>
#include <cando/chem/readAmberParameters.h>
#include <cando/chem/restraint.h>
#include <cando/chem/rotamer.h>
#include <cando/chem/searchStatistics.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/specificContext.h>
#include <cando/chem/statusTracker.h>
#include <cando/chem/stereochemistry.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <cando/chem/structureComparer.h>
#include <cando/chem/structureList.h>
#include <cando/chem/superposeEngine.h>
#include <cando/chem/topology.h>
#include <cando/chem/trajectory.h>
#include <cando/chem/twister.h>
#include <cando/chem/virtualSphere.h>
#include <cando/chem/zMatrix.h>
#include <clasp/core/activationFrame.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/hashTableEql.h>
#include <clasp/core/hashTableEqual.h>
#include <clasp/core/hashTableEqualp.h>
#include <clasp/core/instance.h>
#include <clasp/core/null.h>
#include <clasp/core/singleDispatchGenericFunction.h>
#include <clasp/core/specialForm.h>
#include <cando/geom/boundingBox.h>
#include <cando/geom/color.h>
#include <cando/geom/omatrix.h>
#include <cando/geom/ovector2.h>
#include <cando/geom/ovector3.h>
#include <cando/units/dimension.h>
#include <cando/units/quantity.h>
#include <cando/units/unit.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/atom.h>
#include <cando/chem/constitution.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyFixedNonbond.h>
#include <cando/chem/energyImproperRestraint.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyStretch.h>
#include <cando/chem/entityNameSet.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/ffItorDb.h>
#include <cando/chem/ffPtorDb.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffVdwDb.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/superposableConformationCollection.h>
#include <clasp/core/sexpLoadArchive.h>
#include <clasp/core/sexpSaveArchive.h>
#include <clasp/core/array.h>
#include <cando/geom/coordinateArray.h>
#include <cando/units/namedUnit.h>
#include <cando/chem/virtualAtom.h>
#include <clasp/core/metaClass.h>
#include <cando/chem/monomerPack.h>
#include <cando/chem/representedEntityNameSet.h>
#include <clasp/core/bignum.h>
#include <clasp/clbind/class_rep.h>

namespace adapt
{

};




namespace adapt
{
#if 0
#define EXPOSE_TO_CANDO
#define Use_AdaptPkg
#define EXTERN_REGISTER
#include <clasp/core/initClasses.h>
#undef EXTERN_REGISTER
#undef Use_AdaptPkg
#undef EXPOSE_TO_CANDO
#endif
};





namespace adapt
{

    const char* Adapt_nicknames[] = { "" };


#define SYMBOLS_STATIC
  #ifndef SCRAPING
#include <generated/symbols_scraped_inc.h>
  #endif



    void AdaptExposer_O::expose(core::Lisp_sp lisp,WhatToExpose what) const
    {
	switch (what)
	{
	case candoClasses:
	{
	    /*! Use the core package */
//	    this->package()->usePackage(lisp->findPackage(CorePkg));

	}
	break;
	case candoFunctions:
	{
	    // nothing
	}
	break;
	case candoGlobals:
	{
	}
	break;
	case pythonClasses:
	{
	}
	break;
	case pythonGlobals:
	{
	    // nothing currently
	}
	break;
	case pythonFunctions:
	{
	    // nothing
	}
	break;
	};
    };
    
};





extern "C"
{
    bool ___user_libadapt()
    {
	static adapt::AdaptExposer_O* adaptPkgP = NULL;
	if ( adaptPkgP == NULL )
	{
	    printf("%s:%d - initializing adapt library xxx\n", __FILE__, __LINE__ );
	    adaptPkgP = new adapt::AdaptExposer_O(_lisp);
	    _lisp->installPackage(adaptPkgP);
	    return true;
	} else return false;
    }

}
