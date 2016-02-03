#define	DEBUG_LEVEL_FULL
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



void	Twister_O::initialize()
{
    this->Base::initialize();
    this->_FixedRef = _Nil<Atom_O>();
    this->_Fixed = _Nil<Atom_O>();
    this->_Movable = _Nil<Atom_O>();
    this->_MovableRef = _Nil<Atom_O>();
    this->_Atoms.clear();
}



#ifdef XML_ARCHIVE
    void	Twister_O::archive(core::ArchiveP node)
{
    node->attribute("Fixed",this->_Fixed);
    node->attribute("Movable",this->_Movable);
    node->archiveVector0("MoveAtoms",this->_Atoms);
}
#endif



CL_LISPIFY_NAME("setFixedAndMovable");
CL_DEFMETHOD void	Twister_O::setFixedAndMovable(Atom_sp fixed, Atom_sp movable )
{
    this->_Fixed = fixed;
    this->_Movable = movable;
}

CL_LISPIFY_NAME("addAtom");
CL_DEFMETHOD void	Twister_O::addAtom(Atom_sp a)
{
    this->_Atoms.push_back(a);
}


/*!
 * If force2Mobile is true then the a2 atoms are always the mobile ones
 */
void	Twister_O::_defineForDihedral(Atom_sp a1ref, Atom_sp a1, Atom_sp a2, Atom_sp a2ref, bool force2Mobile )
{
    BondOrder bo = a1->bondOrderTo(a2);
    bool hasBond = a1->isBondedTo(a2);
    if ( hasBond )
    {
	a1->removeBondTo(a2);
    }
		//
		// Accumulate the atoms on each side of the bond
		//
    core::List_sp a1Atoms = _Nil<core::T_O>();
    SpanningLoop_sp span1 = SpanningLoop_O::create();
    span1->setTop(a1);
    while ( span1->advance() )
    {
	a1Atoms = core::Cons_O::create(span1->getAtom(),a1Atoms);
    }
    core::List_sp a2Atoms = _Nil<core::T_O>();
    SpanningLoop_sp span2 = SpanningLoop_O::create();
    span2->setTop(a2);
    while ( span2->advance() )
    {
	a2Atoms = core::Cons_O::create(span2->getAtom(),a2Atoms);
    }
    core::List_sp twistAtoms;
    Atom_sp twistFixed, twistMovable, twistFixedRef, twistMovableRef;
    if ( force2Mobile )
    {
	twistAtoms = a2Atoms;
	twistFixedRef = a1ref;
	twistFixed = a1;
	twistMovable = a2;
	twistMovableRef = a2ref;
    } else
    {
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
    if ( hasBond )
    {
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



CL_LISPIFY_NAME("defineForDihedral");
CL_DEFMETHOD void	Twister_O::defineForDihedral(Atom_sp a1ref, Atom_sp a1, Atom_sp a2, Atom_sp a2ref )
{
    this->_defineForDihedral(a1ref,a1,a2,a2ref,false);
}




CL_LISPIFY_NAME("defineForBond");
CL_DEFMETHOD void	Twister_O::defineForBond(Atom_sp a1, Atom_sp a2 )
{
    this->_defineForDihedral( _Nil<Atom_O>(),a1,a2,_Nil<Atom_O>(),false);
}


CL_LISPIFY_NAME("defineFixedAndMobile");
CL_DEFMETHOD void	Twister_O::defineFixedAndMobile(Atom_sp afixed, Atom_sp amobile )
{
    this->_defineForDihedral(_Nil<Atom_O>(),afixed,amobile,_Nil<Atom_O>(), true);
}


CL_LISPIFY_NAME("rotate");
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




CL_LISPIFY_NAME("rotateAbsolute");
CL_DEFMETHOD void	Twister_O::rotateAbsolute(double angle)
{
    if (this->_FixedRef.nilp() )
    {
	SIMPLE_ERROR(BF("For absolute rotations the reference atoms of the twister must be defined"));
    }
    Vector3 v1 = this->_FixedRef->getPosition();
    Vector3 v2 = this->_Fixed->getPosition();
    Vector3 v3 = this->_Movable->getPosition();
    Vector3 v4 = this->_MovableRef->getPosition();

    double currentAngle = calculateDihedral(v1,v2,v3,v4);
    double delta = angle - currentAngle;
    _lisp->print(BF("current dihedral = %8.3lf  desired dihedral = %8.3lf delta=%8.3lf") % (currentAngle/0.0174533) % (angle/0.0174533) % (delta/0.0174533) );
    this->rotate(delta);
#ifdef DEBUG_ON
    v1 = this->_FixedRef->getPosition();
    v2 = this->_Fixed->getPosition();
    v3 = this->_Movable->getPosition();
    v4 = this->_MovableRef->getPosition();
    double resultAngle = calculateDihedral(v1,v2,v3,v4);
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



#ifdef XML_ARCHIVE
void	TwisterDriver_O::archive(core::ArchiveP node)
{
    node->attribute("Aggregate",this->_Aggregate);
    node->archiveVector0("Twisters",this->_Twisters);
}
#endif


CL_LISPIFY_NAME("addTwister");
CL_DEFMETHOD void	TwisterDriver_O::addTwister(	Twister_sp twister )
{
    this->_Twisters.push_back(twister);
}


CL_LISPIFY_NAME("getTwister");
CL_DEFMETHOD Twister_sp	TwisterDriver_O::getTwister(uint i)
{_OF();
    ASSERT_lessThan(i,this->_Twisters.size());
    return this->_Twisters[i];
}

CL_LISPIFY_NAME("perturbConformation");
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

