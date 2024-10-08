/*
    File: mbbCoreTools.cc
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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define DEBUG_LEVEL_NONE
//
// (C) 2004 Christian E. Schafmeister
//

#define TURN_DEBUG_OFF

#include <clasp/core/common.h>

#include <stdio.h>
#include <string>
#include <math.h>
#include <cando/geom/ovector3.h>
#include <cando/geom/omatrix.h>
#include <clasp/core/exceptions.h>
#include <cando/chem/matter.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/chem/loop.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/mbbCoreTools.h>
#include <clasp/core/wrappers.h>

namespace chem {

#if 0  //[
    Vector3	geometricCenterOfMatter( const Matter_sp& c )
{
    Vector3	sum;
int	count;
Atom_sp	a;
Loop	l;


    sum.set(0,0,0);
    count = 0;
    l.loopTopGoal( c, ATOMS );
    while (l.advanceLoopAndProcess()) {
	a = l.getAtom();
	sum = sum + a->getPosition();
	count++;
    }
    LOG("geometricCenterOfMatter number of atoms = {}" , (count ) );
    if ( count != 0 ) {
        sum = sum*(1.0/count);
    }
    LOG("Geometric center = {}, {}, {}" , (sum.getX()) , (sum.getY()) , (sum.getZ() ) );
    return sum;
}


    Vector3 geometricCenterOfAggregate( const Aggregate_sp& c )
{
    Vector3	v;
    v = geometricCenterOfMatter(c);
    return v;
}
    Vector3 geometricCenterOfMolecule( const Molecule_sp& c )
{
    Vector3	v;
    v = geometricCenterOfMatter(c);
    return v;
}
    Vector3 geometricCenterOfResidue( const Residue_sp& c )
{
    Vector3	v;
    v = geometricCenterOfMatter(c);
    return v;
}
#endif //]

//
//
double closestInterResidueDist2(Residue_sp res1, Residue_sp res2) {
  Loop l;
  gctools::Vec0<Atom_sp> atoms;
  Atom_sp a;
  double dist2, temp;
  double x, y, z;
  gctools::Vec0<Atom_sp>::iterator ai;

  dist2 = std::numeric_limits<double>::max();
  l.loopTopGoal(res2, ATOMS);
  while (l.advanceLoopAndProcess()) {
    a = l.getAtom();
    atoms.push_back(a);
  }
  l.loopTopGoal(res1, ATOMS);
  while (l.advanceLoopAndProcess()) {
    a = l.getAtom();
    for (ai = atoms.begin(); ai != atoms.end(); ai++) {
      x = (*ai)->getPosition().getX() - a->getPosition().getX();
      if (x > dist2)
        continue;
      y = (*ai)->getPosition().getY() - a->getPosition().getY();
      if (y > dist2)
        continue;
      z = (*ai)->getPosition().getZ() - a->getPosition().getZ();
      if (z > dist2)
        continue;
      temp = x * x + y * y + z * z;
      if (temp < dist2)
        dist2 = temp;
    }
  }
  return dist2;
}

//
//	closestNonBondApproachInMolecule
//
//	Return the closest approach between two atoms
//	that are not in contiguous residues
//
double closestNonBondApproachInMolecule(Molecule_sp c) {
  gctools::Vec0<Residue_sp> residues;
  double thisDist2, dist2;
  int rmain, rtemp;
  Loop l;
  Residue_sp res;

  dist2 = std::numeric_limits<double>::max();

  l.loopTopGoal(Matter_sp(c), RESIDUES);
  while ((l.advanceLoopAndProcess())) {
    res = l.getResidue();
    residues.push_back(res);
  }

  // Loop over all the residues
  for (rmain = residues.size() - 1; rmain >= 0; rmain--) {
    if (rmain > 1) {
      for (rtemp = rmain - 2; rtemp >= 0; rtemp--) {
        res = residues[rtemp];
        thisDist2 = closestInterResidueDist2(residues[rmain], res);
        if (thisDist2 < dist2) {
          dist2 = thisDist2;
        }
      }
    }
    if (rmain < ((int)residues.size() - 2)) {
      for (rtemp = rmain + 2; rtemp < (int)residues.size(); rtemp++) {
        res = residues[rtemp];
        thisDist2 = closestInterResidueDist2(residues[rmain], res);
        if (thisDist2 < dist2) {
          dist2 = thisDist2;
        }
      }
    }
  }
  return sqrt(dist2);
}

uint connectVdwOverlappingUnconnectedAtoms(gctools::Vec0<Atom_sp> &atoms) {
  uint connected = 0;
  LOG("There are {} unconnected atoms to connect together", atoms.size());
  for (gctools::Vec0<Atom_sp>::iterator i1 = atoms.begin(); i1 != atoms.end() - 1; i1++) {
    Vector3 v1 = (*i1)->getPosition();
    //	    double rad1 = vdwRadiusForElement((*i1)->getElement());
    for (gctools::Vec0<Atom_sp>::iterator i2 = i1 + 1; i2 != atoms.end(); i2++) {
      Vector3 v2 = (*i2)->getPosition();
      //		double rad2 = vdwRadiusForElement((*i2)->getElement());
      Vector3 delta = v2.sub(v1);
      double len = delta.length();
      double maxLen = maximumBondLengthBetweenElements((*i1)->getElement(), (*i2)->getElement());
      if (len < maxLen) {
        if (!(*i1)->isBondedTo(*i2)) {
          LOG("Forming bond between atom({}) and atom({}) - they were {} angstroms apart < max({})", (*i1)->description(),
              (*i2)->description(), len, maxLen);
          (*i1)->bondToSingle(*i2);
          connected++;
        }
      }
    }
  }
  LOG("Connected {} bonds", connected);
  return connected;
}

DOCGROUP(cando);
CL_DEFUN void chem__connectAtomsInMatterInCovalentContact(Matter_sp matter) {
  gctools::Vec0<Atom_sp> atoms;
  Loop lAtoms;
  lAtoms.loopTopGoal(matter, ATOMS);
  while (lAtoms.advance()) {
    atoms.push_back(lAtoms.getAtom());
  }
#ifdef DEBUG_ON
  uint numVdwConnected = connectVdwOverlappingUnconnectedAtoms(atoms);
  LOG("Number of atoms connected using vdw overlap: {}", numVdwConnected);
#endif
}

/*! Return all of the atoms in the matter that are within the sphere indicated by center/radius
 */
DOCGROUP(cando);
CL_DEFUN adapt::ObjectSet_sp chem__atomsWithinSphereAsObjectSet(Matter_sp matter, Vector3 center, double radius) {
  LOG("Looking for atoms in sphere with center({}) and radius({})", center.asString(), radius);
  adapt::ObjectSet_sp result = adapt::ObjectSet_O::create();
  Loop lAtoms;
  lAtoms.loopTopGoal(matter, ATOMS);
  while (lAtoms.advance()) {
    Atom_sp a = lAtoms.getAtom();
    Vector3 pos = a->getPosition();
    Vector3 delta = pos.sub(center);
    double dist = delta.length();
    LOG("Looking at atom at position({}) distance({}) insphere({})", pos.asString(), dist, (dist < radius));
    if (dist < radius) {
      result->insert(a);
    }
  }
  return result;
}

/*! Return all of the atoms in the matter that are within the sphere indicated by center/radius
 */
DOCGROUP(cando);
CL_DEFUN adapt::ObjectSet_sp chem__atomsInMatterAsObjectSet(Matter_sp matter) {
  adapt::ObjectSet_sp result = adapt::ObjectSet_O::create();
  Loop lAtoms;
  lAtoms.loopTopGoal(matter, ATOMS);
  while (lAtoms.advance()) {
    Atom_sp a = lAtoms.getAtom();
    result->insert(a);
  }
  return result;
}

}; // namespace chem
