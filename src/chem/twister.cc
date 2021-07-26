/*
    File: twister.cc
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
#define	DEBUG_LEVEL_NONE
//
// (C) 2004 Christian E. Schafmeister
//



#include <clasp/core/common.h>

#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <exception>
#include <cando/chem/matter.h>
#include <clasp/core/exceptions.h>
#include <cando/chem/loop.h>
#include <cando/chem/spanningLoop.h>
//#include	"xmlGraphics.h"
#include <cando/chem/twister.h>
#include <clasp/core/numerics.h>
#include <cando/chem/spanningLoop.h>
#include <clasp/core/wrappers.h>

namespace chem {


CL_DOCSTRING(R"doc(Create a new twister instance.)doc");
CL_DEFUN Twister_sp chem__make_twister()
{
  auto td = gctools::GC<Twister_O>::allocate_with_default_constructor();
  return td;
}

void	Twister_O::initialize()
{
    this->Base::initialize();
    this->_FixedRef = _Nil<Atom_O>();
    this->_Fixed = _Nil<Atom_O>();
    this->_Movable = _Nil<Atom_O>();
    this->_MovableRef = _Nil<Atom_O>();
    this->_Atoms.clear();
}


CL_LAMBDA(twister fixed-atom movable-atom);
CL_DOCSTRING(R"doc(Set the fixed and moveable atoms of the twister.)doc");
CL_LISPIFY_NAME("Twister_setFixedAndMovable");
CL_DEFMETHOD void	Twister_O::setFixedAndMovable(Atom_sp fixed, Atom_sp movable )
{
    this->_Fixed = fixed;
    this->_Movable = movable;
}

CL_LAMBDA(twister atom);
CL_DOCSTRING(R"doc(Add an atom to the twister.)doc");
CL_LISPIFY_NAME("Twister_addAtom");
CL_DEFMETHOD void	Twister_O::addAtom(Atom_sp a)
{
    this->_Atoms.push_back(a);
}


/*!
 * If force2Mobile is true then the a2 atoms are always the mobile ones.
 * Otherwise use the side of the bond that has the fewest atoms as the mobile ones.
 */
void	Twister_O::_defineForDihedral(Atom_sp a1ref, Atom_sp a1, Atom_sp a2, Atom_sp a2ref, bool force2Mobile )
{
    BondOrder bo = a1->bondOrderTo(a2);
    bool hasBond = a1->isBondedTo(a2);
    if ( hasBond ) {
      a1->removeBondTo(a2);
    }
		//
		// Accumulate the atoms on each side of the bond
		//
    core::List_sp a1Atoms = _Nil<core::T_O>();
    SpanningLoop_sp span1 = SpanningLoop_O::create();
    span1->setTop(a1);
    while ( span1->advance() ) {
	a1Atoms = core::Cons_O::create(span1->getAtom(),a1Atoms);
    }
    core::List_sp a2Atoms = _Nil<core::T_O>();
    SpanningLoop_sp span2 = SpanningLoop_O::create();
    span2->setTop(a2);
    while ( span2->advance() ) {
	a2Atoms = core::Cons_O::create(span2->getAtom(),a2Atoms);
    }
    if (a1ref.nilp()) {
      // Find the first heavy atom connected to a1 and make it a1ref
      for (int ii=0; ii<a1->numberOfBonds(); ii++ ) {
        Atom_sp other = a1->bondedNeighbor(ii);
        if (other->getAtomicNumber()>1) {
          a1ref = other;
          break;
        }
      }
    }
    if (a2ref.nilp()) {
      // Find the first heavy atom connected to a2 and make it a2ref
      for (int ii=0; ii<a2->numberOfBonds(); ii++ ) {
        Atom_sp other = a2->bondedNeighbor(ii);
        if (other->getAtomicNumber()>1) {
          a2ref = other;
          break;
        }
      }
    }
    core::List_sp twistAtoms;
    Atom_sp twistFixed, twistMovable, twistFixedRef, twistMovableRef;
    if ( force2Mobile ) {
	twistAtoms = a2Atoms;
	twistFixedRef = a1ref;
	twistFixed = a1;
	twistMovable = a2;
	twistMovableRef = a2ref;
    } else {
      if ( core::cl__length(a1Atoms) < core::cl__length(a2Atoms) )
	{
	    twistAtoms = a1Atoms;
	    twistFixedRef = a2ref;
	    twistFixed = a2;
	    twistMovable = a1;
	    twistMovableRef = a1ref;
	} else
	{
	    twistAtoms = a2Atoms;
	    twistFixedRef = a1ref;
	    twistFixed = a1;
	    twistMovable = a2;
	    twistMovableRef = a2ref;
	}
    }
    	//
	// Restore bond
	//
    if ( hasBond ) {
	a1->bondTo(a2,bo);
    }
    this->_Atoms.clear();
    this->setFixedAndMovable(twistFixed,twistMovable);
    this->_FixedRef = twistFixedRef;
    this->_MovableRef = twistMovableRef;
    for ( auto cur : twistAtoms ) {
	this->addAtom(cur->car<Atom_O>());
    }
}



CL_LAMBDA(twister atom1-ref atom1 atom2 atom2-ref);
CL_DOCSTRING(R"doc(Define the twister based on a dihedral between four atoms. Whichever side of the
dihedral contains the least atoms will be the movable side while the other side
will remain fixed.)doc");
CL_LISPIFY_NAME("twister_defineForDihedral");
CL_DEFMETHOD void	Twister_O::defineForDihedral(Atom_sp a1ref, Atom_sp a1, Atom_sp a2, Atom_sp a2ref )
{
    this->_defineForDihedral(a1ref,a1,a2,a2ref,false);
}




CL_LAMBDA(twister atom1 atom2);
CL_DOCSTRING(R"doc(Define the twister based on a bond. Whichever side of the bond contains the
least atoms will be the movable side while the other side will remain fixed.)doc");
CL_LISPIFY_NAME("twister_defineForBond");
CL_DEFMETHOD void	Twister_O::defineForBond(Atom_sp a1, Atom_sp a2 )
{
    this->_defineForDihedral( _Nil<Atom_O>(),a1,a2,_Nil<Atom_O>(),false);
}


CL_LAMBDA(twister fixed-atom movable-atom);
CL_DOCSTRING(R"doc(Define the twister based on a fixed and movable atom.)doc");
CL_LISPIFY_NAME("twister_defineFixedAndMobile");
CL_DEFMETHOD void	Twister_O::defineFixedAndMobile(Atom_sp afixed, Atom_sp amobile )
{
    this->_defineForDihedral(_Nil<Atom_O>(),afixed,amobile,_Nil<Atom_O>(), true);
}


CL_DOCSTRING(R"doc(Rotate the twister by a relative angle in radians.)doc");
CL_LAMBDA(twister angle);
CL_LISPIFY_NAME("twister_rotate");
CL_DEFMETHOD void	Twister_O::rotate(double angle)
{_OF();
    Matrix	transform, tm;
    Vector3	rotVec, tv;

    LOG(BF("Rotating by %lf degrees") % (angle/0.0174533) );
    ASSERTNOTNULLP(this->_Fixed,"Fixed atom undefined");
    ASSERTNOTNULLP(this->_Movable,"Movable atom is undefined");
    ASSERTP(this->_Atoms.size()>0,"There must be atoms to rotate");
    rotVec = this->_Movable->getPosition()-this->_Fixed->getPosition();
    rotVec = rotVec.normalized();
    tv = this->_Movable->getPosition();
    tv = tv.multiplyByScalar(-1.0);
    transform.translate(tv);
    LOG(BF("To origin transform matrix: %s") % (transform.asString().c_str()) );
    tm.rotationAxis(angle,&rotVec);
    transform = tm*transform;
    tv = this->_Movable->getPosition();
    tm.translate(tv);
    LOG(BF("From origin transform matrix: %s") % (tm.asString().c_str()) );
    transform = tm*transform;
    LOG(BF("Complete transform matrix: %s") % (transform.asString().c_str()) );
    for ( gctools::Vec0<Atom_sp>::iterator ai=this->_Atoms.begin();
	  ai!=this->_Atoms.end(); ai++ ) {
	LOG(BF("Perturbing atom(%s)  start position = %lf, %lf, %lf") % (*ai)->description().c_str() % (*ai)->getPosition().getX() % (*ai)->getPosition().getY() % (*ai)->getPosition().getZ()  );
	tv = (*ai)->getPosition();
	tv = transform*tv;
	(*ai)->setPosition(tv);
	LOG(BF("                     pert. position = %lf, %lf, %lf") % (*ai)->getPosition().getX() % (*ai)->getPosition().getY() % (*ai)->getPosition().getZ()  );
    }
}


CL_LAMBDA(twister);
CL_DOCSTRING(R"doc(Get the current dihedral angle of the twister in radians.)doc");
CL_LISPIFY_NAME("twister_currentDihedralAngleRadians");
CL_DEFMETHOD double Twister_O::currentDihedralAngleRadians()
{
  if (this->_FixedRef.nilp() )
  {
    SIMPLE_ERROR(BF("For absolute rotations the reference atoms of the twister must be defined"));
  }
  Vector3 v1 = this->_FixedRef->getPosition();
  Vector3 v2 = this->_Fixed->getPosition();
  Vector3 v3 = this->_Movable->getPosition();
  Vector3 v4 = this->_MovableRef->getPosition();
  
  double currentAngle = geom::calculateDihedral(v1,v2,v3,v4);
  return currentAngle;
}



CL_LAMBDA(twister angle);
CL_DOCSTRING(R"doc(Rotate the twister to an absolute angle in radians.)doc");
CL_LISPIFY_NAME("twister_rotateAbsolute");
CL_DEFMETHOD void Twister_O::rotateAbsolute(double angle)
{
  double currentAngle = this->currentDihedralAngleRadians();
  double delta = angle - currentAngle;
  if (chem__verbose(1)) {
    _lisp->print(BF("current dihedral = %8.3lf  desired dihedral = %8.3lf delta=%8.3lf") % (currentAngle/0.0174533) % (angle/0.0174533) % (delta/0.0174533) );
  }
  this->rotate(delta);
#ifdef DEBUG_ON
  v1 = this->_FixedRef->getPosition();
  v2 = this->_Fixed->getPosition();
  v3 = this->_Movable->getPosition();
  v4 = this->_MovableRef->getPosition();
  double resultAngle = geom::calculateDihedral(v1,v2,v3,v4);
  _lisp->print(BF("      result dihedral = %8.3lf") % (resultAngle/0.0174533) );
#endif
}


string	Twister_O::description() const
{
stringstream	ss;
    ss << "Twister[Fixed[" << this->_Fixed->description();
    ss << "] Moveable[" << this->_Movable->description() << "]";
    return ss.str();
}


CL_LAMBDA(aggregate);
CL_DOCSTRING(R"doc(Create a new twister driver based on an aggregate.)doc");
CL_DEFUN TwisterDriver_sp chem__make_twister_driver(Aggregate_sp agg)
{
  auto td = gctools::GC<TwisterDriver_O>::allocate(agg);
  return td;
}


CL_LAMBDA(twister-driver twister);
CL_DOCSTRING(R"doc(Add twister to a twister driver.)doc");
CL_LISPIFY_NAME("twister_driver_addTwister");
CL_DEFMETHOD void	TwisterDriver_O::addTwister(	Twister_sp twister )
{
    this->_Twisters.push_back(twister);
}


CL_LAMBDA(twister-driver index);
CL_DOCSTRING(R"doc(Get the twister in a twister driver based on index.)doc");
CL_LISPIFY_NAME("twister_driver_getTwister");
CL_DEFMETHOD Twister_sp	TwisterDriver_O::getTwister(uint i)
{_OF();
    ASSERT_lessThan(i,this->_Twisters.size());
    if (i>=0 && i<this->_Twisters.size()) return this->_Twisters[i];
    SIMPLE_ERROR(BF("Illegal twister index %d there are only %d twisters") % i % this->_Twisters.size());
}

CL_LAMBDA(twister-driver);
CL_DOCSTRING(R"doc(Randomly select a twister and perturb the angle by a random amount.)doc");
CL_LISPIFY_NAME("twister_driver_perturbConformation");
CL_DEFMETHOD void	TwisterDriver_O::perturbConformation()
{_OF();
int	it,ia;
double	ang;
    if ( this->_Twisters.size() == 0 )
    {
        LOG(BF("There are no twisters and nothing to perturb") );
	return;
    }
    it = core::randomNumber01()*this->_Twisters.size();
    ia = (int)(core::randomNumber01()*12.0);
    ang = (30.0*ia)*0.0174533+(core::randomNumberNormal01()*15);
    LOG(BF("Perturbing twister#%d of %d by %lf degrees") % it % this->_Twisters.size() % (ang/0.0174533)  );
    this->_Twisters[it]->rotate(ang);
}

string	TwisterDriver_O::description() const
{
stringstream	ss;
    ss << "TwisterDriver[#twisters("<<this->_Twisters.size()<<")]";
    return ss.str();
}

};
