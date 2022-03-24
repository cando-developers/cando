/*
    File: plug.cc
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

#include <clasp/core/common.h>
#include <clasp/core/record.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/plug.h>
#include <cando/chem/candoDatabase.h>
#include <clasp/core/cons.h>
#include <cando/chem/coupling.h>
#include <clasp/core/environment.h>
#include <clasp/core/numbers.h>
#include <cando/chem/constitution.h>
#include <clasp/core/wrappers.h>


namespace chem
{
core::Symbol_sp RingClosingMate_O::getName() const
{
  return _lisp->intern("UnnamedRingClosingMate");
};

void Mate_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,cap),this->_Cap);
  this->Base::fields(node);
}


/*
  __BEGIN_DOC(classes.RingClosingMate)
  \inheritsFrom{EntityNameSetWithCap}
  A EntityNameSet that keeps track of a capping monomer that is used to cap training oligomers when they are being defined.
  The capping monomer is supposed to be small and best represent the other members of the RingClosingMate.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.RingClosingMate.!class.RingClosingMate)
  \requiredKeyed{groupNames:}{List::groupNames}

  Initialize a RingClosingMate object.
  __END_DOC
*/

core::Symbol_sp Mate_O::getName() const
{
  return nil<core::Symbol_O>();
}


/*
  __BEGIN_DOC(classes.Mate)
  \inheritsFrom{EntityNameSetWithCap}
  A EntityNameSet that keeps track of a capping monomer that is used to cap training oligomers when they are being defined.
  The capping monomer is supposed to be small and best represent the other members of the Mate.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.Mate.!class.Mate)
  \requiredKeyed{cap:}{Text::capName}
  \requiredKeyed{groupNames:}{List::groupNames}

  Initialize a Mate object.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.Plug)
  Defines one or two atoms of this monomer that can be plugged into, a plug name and a
  collection of Mate objects that can act as mates for this plug.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.Plug.!class.Plug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a Plug object. Plugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  Outgoing plugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  __END_DOC
*/

void Plug_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,name),this->_Name);
}


string Plug_O::__repr__() const
{
  stringstream ss;
  ss << "#<" << this->className() << " " << _rep_(this->_Name) << ">";
  return ss.str();
}

CL_NAME(chem:|Plug_O::getB0|);
CL_DEFMETHOD core::Symbol_sp Plug_O::getB0() const	{ SUBIMP(); };
CL_NAME(chem:|Plug_O::getB1|);
CL_DEFMETHOD core::Symbol_sp Plug_O::getB1() const	{ SUBIMP(); };
CL_NAME(chem:|Plug_O::getBondOrder0|);
CL_DEFMETHOD BondOrder Plug_O::getBondOrder0() const { SUBIMP(); };
CL_NAME(chem:|Plug_O::getBondOrder1|);
CL_DEFMETHOD BondOrder Plug_O::getBondOrder1() const { SUBIMP(); };

CL_NAME(chem:|Plug_O::setName|);
CL_DEFMETHOD void Plug_O::setName(core::Symbol_sp s) {
  this->_Name = s;
}

CL_NAME(chem:|Plug_O::getName|);
CL_DEFMETHOD core::Symbol_sp Plug_O::getName() const {
  return this->_Name;
}

CL_DEFMETHOD bool Plug_O::getIsIn() const {_OF(); SUBCLASS_MUST_IMPLEMENT();};
CL_DEFMETHOD bool Plug_O::getIsOrigin() const {return false;};
CL_DEFMETHOD bool Plug_O::getIsRingClosing() const { return false; };



CL_DEFMETHOD bool Plug_O::getIsOut() const
{
  return !(this->getIsIn()||this->getIsRingClosing());
}

CL_LISPIFY_NAME("otherSidePlugName");
CL_DEFMETHOD     core::Symbol_sp Plug_O::otherSidePlugName()
{
  return DirectionalCoupling_O::otherPlugName(this->_Name);
}


/*
  __BEGIN_DOC(classes.PlugWithMates)
  Defines one or two atoms of this monomer that can be plugged into, a plug name and a
  collection of Mate objects that can act as mates for this plug.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.PlugWithMates.!class.PlugWithMates)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a PlugWithMates object. PlugWithMatess can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  __END_DOC
*/

void PlugWithMates_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,b0),this->_B0);
  node->field_if_not_nil(INTERN_(kw,b1),this->_B1);
  node->field(INTERN_(kw,mates),this->_Mates);
  this->Base::fields(node);
}

  
bool	PlugWithMates_O::recognizesMateNameOrPdb(core::Symbol_sp name)
{
  gctools::Vec0<Mate_sp>::iterator	mi;
  for ( mi=this->_Mates.begin(); mi!=this->_Mates.end(); mi++ ) {
    if ( (*mi)->recognizesNameOrPdb(name) ) return true;
  }
  return false;
}


/*
  __BEGIN_DOC(classes.OutPlug.!class.OutPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a OutPlug object. OutPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  Outgoing plugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  __END_DOC
*/

void OutPlug_O::fields(core::Record_sp node)
{
  node->field_if_not_nil(INTERN_(kw,stubPivotAtom),this->_StubPivotAtom);
  node->field_if_not_default(INTERN_(kw,ringClosing),this->_IsRingClosing,false);
  this->Base::fields(node);
}

CL_LISPIFY_NAME("ringClosingMatesAsList");
CL_DEFMETHOD     core::List_sp OutPlug_O::ringClosingMatesAsList()
{_OF();
  return core::Cons_O::createFromVec0(this->_Mates);
}

bool OutPlug_O::recognizesRingClosingMate(core::Symbol_sp mateName)
{_OF();
  gctools::Vec0<Mate_sp>::iterator it;
  for ( it=this->_Mates.begin(); it!=this->_Mates.end(); it++ )
  {
    if ( (*it)->recognizesMonomerName(mateName)) return true;
  }
  return false;
}

/*
  __BEGIN_DOC(classes.InPlug.!class.InPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a InPlug object. InPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  __END_DOC
*/
void InPlug_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
}


/*
  __BEGIN_DOC(classes.JumpPlug.!class.JumpPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{originFrame:}{Frame::originFrame}

  Initialize a JumpPlug object. JumpPlugs don't make bonds don't have mates but they do have an origin frame that is within the topology.
  __END_DOC
*/

void JumpPlug_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,jumpPointName),this->_JumpPointName);
  this->Base::fields(node);
}


void OriginPlug_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,originAtomName),this->_OriginAtomName);
  this->Base::fields(node);
}

/*
  __BEGIN_DOC(classes.RingClosingPlug.!class.RingClosingPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}
  \optionalKeyed{exportFrame:}{Frame::exportFrame}
  \requiredKeyed{ringClosingMates:}{Cons::listOfRingClosingMates}

  Initialize a RingClosingPlug object. RingClosingPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  RingClosingPlugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  There is an additional list of ringClosingMates that can be attached to this plug without building the mates.
  __END_DOC
*/

#if 0
void RingClosingPlug_O::fields(core::Record_sp node) {
  this->Base::fields(node);
}
#endif

#if 0
RingClosingPlug_O::RingClosingPlug_O(const RingClosingPlug_O& p) : RingClosingPlug_O::Base(p)
{
  LOG(BF("Original %s") % p.description().c_str()  );
  gctools::Vec0<RingClosingMate_sp>::const_iterator vi;
  LOG(BF("Copying %d ring closing mates") % p._Mates.size()  );
  for ( vi=p._RingClosingMates.begin(); vi!=p._RingClosingMates.end(); vi++ ) {
    LOG(BF("Copied mate") );
    auto  rn = gctools::GC<RingClosingMate_O>::copy( *vi->get() ); // = RP_Copy<RingClosingMate_O>(*vi);
    this->_RingClosingMates.push_back(rn);
  }
  LOG(BF("Copy constructed %s") % this->description().c_str()  );
}
#endif

#if 0

CL_LISPIFY_NAME("ringClosingMatesAsList");
CL_DEFMETHOD     core::List_sp RingClosingPlug_O::ringClosingMatesAsList()
{_OF();
  return core::Cons_O::createFromVec0(this->_Mates);
}


bool RingClosingPlug_O::recognizesRingClosingMate(core::Symbol_sp mateName)
{_OF();
  gctools::Vec0<Mate_sp>::iterator it;
  for ( it=this->_Mates.begin(); it!=this->_Mates.end(); it++ )
  {
    if ( (*it)->recognizesMonomerName(mateName)) return true;
  }
  return false;
}

#endif


}; //namespace chem
