       
       
//
// (C) 2004 Christian E. Schafmeister
//



#ifndef	mbbCoreTools_H
#define mbbCoreTools_H

#include <string>
#include <clasp/core/common.h>
#include <cando/chem/matter.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/adapt/objectSet.h>
#include <cando/geom/vector3.h>
#include <cando/geom/matrix.h>

#include <cando/chem/chemPackage.h>


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
