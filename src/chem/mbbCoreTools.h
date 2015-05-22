       
       
//
// (C) 2004 Christian E. Schafmeister
//



#ifndef	mbbCoreTools_H
#define mbbCoreTools_H

#include <string>
#include "core/common.h"
#include "matter.h"
#include "aggregate.h"
#include "molecule.h"
#include "residue.h"
#include "atom.h"
#include "addon/vector3.h"
#include "addon/matrix.h"

#include "chemPackage.h"

#include "core/objectSet.fwd.h"// mbbCoreTools.h wants ObjectSet needs objectSet.fwd.h

namespace chem {


double	closestNonBondApproachInMolecule( Molecule_sp c );
double	closestInterResidueDist2( Residue_sp res1, Residue_sp res2 );


void	testStuff();

void connectAtomsInMatterInCovalentContact(Matter_sp matter);

SMART(ObjectSet);

core::ObjectSet_sp atomsWithinSphere(Matter_sp matter, Vector3 center, double radius);


core::ObjectSet_sp atomsInMatter(Matter_sp matter);


void exposeMbbCoreTools(core::Lisp_sp lisp);

};

#endif
