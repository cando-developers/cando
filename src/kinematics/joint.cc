/*
    File: joint.cc
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
#define DEBUG_LEVEL_FULL


#include <clasp/core/foundation.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/record.h>
#include <clasp/core/ql.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/joint.h>


SYMBOL_EXPORT_SC_(KeywordPkg,noop_to_external);
SYMBOL_EXPORT_SC_(KeywordPkg,noop_to_internal);
SYMBOL_EXPORT_SC_(KeywordPkg,origin_to_external);
SYMBOL_EXPORT_SC_(KeywordPkg,origin_to_internal);
SYMBOL_EXPORT_SC_(KeywordPkg,general_to_external);
SYMBOL_EXPORT_SC_(KeywordPkg,general_to_internal);
SYMBOL_EXPORT_SC_(KeywordPkg,bond_to_external);
SYMBOL_EXPORT_SC_(KeywordPkg,bond_to_internal);
SYMBOL_EXPORT_SC_(KeywordPkg,bond_angle_to_external);
SYMBOL_EXPORT_SC_(KeywordPkg,bond_angle_to_internal);
SYMBOL_EXPORT_SC_(KeywordPkg,bond_angle_dihedral_to_external);
SYMBOL_EXPORT_SC_(KeywordPkg,bond_angle_dihedral_to_internal);

namespace translate
{
template <>
struct	to_object<kinematics::CoordinateCalculator>
{
  typedef	core::Symbol_sp	ExpectedType;
  typedef	core::Symbol_sp	DeclareType;
  static core::T_sp convert(kinematics::CoordinateCalculator v)
  {_G();
    if (v== kinematics::noop_to_external) {
      return kw::_sym_noop_to_external;
    } else if (v== kinematics::noop_to_internal) {
      return kw::_sym_noop_to_internal;
    } else if (v== kinematics:: origin_to_external) {
      return kw::_sym_origin_to_external;
    } else if (v== kinematics::origin_to_internal) {
      return kw::_sym_origin_to_internal;
    } else if (v== kinematics::general_to_external) {
      return kw::_sym_general_to_external;
    } else if (v== kinematics::general_to_internal) {
      return kw::_sym_general_to_internal;
    } else if (v== kinematics::bond_to_external) {
      return kw::_sym_bond_to_external;
    } else if (v== kinematics::bond_to_internal) {
      return kw::_sym_bond_to_internal;
    } else if (v== kinematics::bond_angle_to_external) {
      return kw::_sym_bond_angle_to_external;
    } else if (v== kinematics::bond_angle_to_internal) {
      return kw::_sym_bond_angle_to_internal;
    } else if (v== kinematics::bond_angle_dihedral_to_external) {
      return kw::_sym_bond_angle_dihedral_to_external;
    } else if (v== kinematics::bond_angle_dihedral_to_internal) {
      return kw::_sym_bond_angle_dihedral_to_internal;
    }
    SIMPLE_ERROR("Cannot convert CoordinateCalculator");
  };
};
template <>
struct translate::from_object<kinematics::CoordinateCalculator>
{
  typedef	kinematics::CoordinateCalculator ExpectedType;
  typedef	ExpectedType DeclareType;
  DeclareType _v;
  from_object(core::T_sp o)
  {
    if (o == kw::_sym_noop_to_external) {
      this->_v = kinematics::noop_to_external;
    } else if (o == kw::_sym_noop_to_internal) {
      this->_v = kinematics::noop_to_internal;
    } else if (o == kw::_sym_origin_to_external) {
      this->_v = kinematics::origin_to_external;
    } else if (o == kw::_sym_origin_to_internal) {
      this->_v = kinematics::origin_to_internal;
    } else if (o == kw::_sym_general_to_external) {
      this->_v = kinematics::general_to_external;
    } else if (o == kw::_sym_general_to_internal) {
      this->_v = kinematics::general_to_internal;
    } else if (o == kw::_sym_bond_to_external) {
      this->_v = kinematics::bond_to_external;
    } else if (o == kw::_sym_bond_to_internal) {
      this->_v = kinematics::bond_to_internal;
    } else if (o == kw::_sym_bond_angle_to_external) {
      this->_v = kinematics::bond_angle_to_external;
    } else if (o == kw::_sym_bond_angle_to_internal) {
      this->_v = kinematics::bond_angle_to_internal;
    } else if (o == kw::_sym_bond_angle_dihedral_to_external) {
      this->_v = kinematics::bond_angle_dihedral_to_external;
    } else if (o == kw::_sym_bond_angle_dihedral_to_internal) {
      this->_v = kinematics::bond_angle_dihedral_to_internal;
    } else {
      SIMPLE_ERROR("Cannot convert CoordinateCalculator {}" , _rep_(o));
    }
  }
};

  

};



namespace kinematics {
FORWARD(JumpJoint);

#define	ASSERT_VALID_HANDLE(tree,handle)				\
  ASSERTF((int)handle<tree->numberOfEntries(),("The handle[{}] is out of range (0->{}]") , handle , tree->numberOfEntries()); \
  ASSERTF(tree->_AtomHolders[handle]._Type != unused,("The handle represents an unused node"));


string Joint_O::__repr__() const {
  stringstream ss;
  ss << "#<";
  ss << this->className();
  ss << " " << _rep_(this->_Name) << ">";
  return ss.str();
}
 
void Joint_O::fields(core::Record_sp node) {
  node->field_if_not_unbound(INTERN_(kw,parent),this->_Parent);
  node->field(INTERN_(kw,name),this->_Name); // name
  node->field(INTERN_(kw,id),this->_Id);
  node->field(INTERN_(kw,posIndexX3),this->_PositionIndexX3);
  node->field(INTERN_(kw,endPosIndexX3),this->_EndPositionIndexX3);
  node->field(INTERN_(kw,properties),this->_Properties);
}

CL_NAME(KIN:JOINT/NAME);
CL_DEFMETHOD core::T_sp Joint_O::name() const {
  return this->_Name;
}


string Joint_O::asString() const
{
  stringstream ss;
  ss << "Node[" << _rep_(this->typeSymbol()) << "]  id=" << this->_Id.asString() << std::endl;
  ss << "Children: ";
  if ( this->_numberOfChildren() == 0 )
  {
    ss << "-none-";
  } else {
    for ( int i=0; i<(int)this->_numberOfChildren(); i++ )
    {
      Joint_sp atom = this->_child(i);
      ss << _rep_(atom);
    }
  }
  ss << std::endl;
  ss << "PositionIndexX3: " << this->_PositionIndexX3 << std::endl;
  return ss.str();
}

void Joint_O::setParent(Joint_sp parent)
{
  this->_Parent = parent;
}
    
    

CL_DEFMETHOD
void Joint_O::insertChild(int before, Joint_sp child )
{
  this->_insertChild(before,child);
  child->setParent(this->asSmartPtr());
}

CL_DEFMETHOD
void Joint_O::appendChild(Joint_sp child)
{
  if ( gc::IsA<JumpJoint_sp>(child) )
  {
    int idx = this->firstNonJumpChildIndex();
    this->_insertChild(idx,child);
    child->setParent(this->asSmartPtr());
    return;
  }
  this->_appendChild(child);
  child->setParent(this->asSmartPtr());
}

CL_DEFMETHOD
void Joint_O::eraseChild(Joint_sp child)
{
  Joint_sp atom = this->asSmartPtr();
  for (int i=0; i<atom->_numberOfChildren(); i++ ) {
    if (child == atom->_child(i)) {
      atom->_releaseChild(i);
      return;
    }
  }
  SIMPLE_ERROR("Could not find child {} in {} - so it could not be erased" , _rep_(child) , _rep_(this->asSmartPtr()));
}


CL_DEFMETHOD
int Joint_O::indexOfChild(Joint_sp child)
{
  Joint_sp atom = this->asSmartPtr();
  for (int i=0; i<atom->_numberOfChildren(); i++ )
  {
    if (child == atom->_child(i)) return i;
  }
  SIMPLE_ERROR("Could not find child");
}



CL_NAME(KIN:JOINT/ADD-CHILD);
CL_DEFMETHOD
void Joint_O::addChild(Joint_sp child)
{
  LOG(("Inserting child: %s") , _rep_(child));
  if ( gc::IsA<JumpJoint_sp>(child) ) {
    LOG(("It's a jump, inserting it at the start"));
    this->_insertChild(0,child);
    child->setParent(this->asSmartPtr());
  } else {
    LOG("It's a non-jump atom");
    LOG("We are at the end of the Children - appending");
    this->_appendChild(child);
    child->setParent(this->asSmartPtr());
  }
}

RootJointInfo const* Joint_O::rootJointInfo() const
{
  THROW_HARD_ERROR(("rootJointInfo is not available for non-root atoms"));
}

int Joint_O::firstNonJumpChildIndex() const
{
  for ( int i=0; i<this->_numberOfChildren(); i++ )
  {
    if ( !gc::IsA<JumpJoint_sp>(this->_child(i)) ) return i;
  }
  return this->_numberOfChildren();
}




void Joint_O::recursiveDumpChildrenIntoStringStream(const string& prefix, stringstream& out)
{
  out << prefix << this->typeSymbol() << "/";
  out << this->_Id.asString() << " ";
#if DEBUG_KIN_ATOM
  out << this->_Comment;
#endif
  out << std::endl;
  string newPrefix = prefix + " ";
  for ( int i=0; i<this->_numberOfChildren(); i++ )
  {
    Joint_sp child = this->_child(i);
    child->recursiveDumpChildrenIntoStringStream(newPrefix,out);
  }
}


bool Joint_O::stubDefined() const
{
	// Inspired by rosetta::Joint_::stub_defined()
	// We have to handle a couple of cases
	//
	// note -- in counting dependent atoms, exclude JumpJoints
	//
	// 1. no dependent atoms -> no way to define new coord sys
	//    on this end. So take parent's rotation and my position
	//
	// 2. one dependent atom -> no way to define unique coordinate system
	//    on this end. So take parents rotation and my position
	//
	// 3. two or more dependent atoms
	//	a) if my first atom has a dependent atom, 
	//	   use myself, my first atom and his first atom
	//	b) otherwise, use myself, my first atom, my second atom
	//
  if ( gc::IsA<JumpJoint_sp>(this->asSmartPtr()))
  {
    Joint_sp first = this->getNonJumpJoint(0);
    if ( first.boundp() &&
         (first->getNonJumpJoint(0) || this->getNonJumpJoint(1).boundp() ))
    {
      return true;
    } else
    {
      return false;
    }
  }
  return true;
}


Joint_sp Joint_O::getNonJumpJoint(int offset) const
{
  int idx = this->firstNonJumpChildIndex();
  idx += offset;
  if ( idx >= this->_numberOfChildren() ) return unbound<Joint_O>();
  return this->_child(idx);
}

CL_DEFMETHOD core::List_sp Joint_O::jointChildren() const
{
  ql::list l;
  for ( int i(0), iEnd(this->_numberOfChildren()); i<iEnd; ++i ) {
    l << this->_child(i);
  }
  return l.cons();
}

CL_DEFUN size_t kin__numberOfChildren(Joint_sp joint) {
  return joint->_numberOfChildren();
}

CL_DEFMETHOD Joint_sp Joint_O::onlyOtherChild(Joint_sp child) const
{
  if (this->_numberOfChildren()!=2) {
    SIMPLE_ERROR("There must be two children but there are {}", this->_numberOfChildren());
  }
  if (this->_child(0) == child) return this->_child(1);
  if (this->_child(1) == child) return this->_child(0);
  SIMPLE_ERROR("Could not find child {} in the two children ({} and {})",
               _rep_(child), _rep_(this->_child(0)), _rep_(this->_child(1)));
}

CL_DEFMETHOD Joint_sp Joint_O::childWithName(core::T_sp name) const
{
  for ( size_t ii=0; ii<this->_numberOfChildren(); ii++ ) {
    if (this->_child(ii)->name() == name) return this->_child(ii);
  }
  SIMPLE_ERROR("Could not find child with name {}", _rep_(name));
}

    /*! See file:///Users/meister/Development/rosetta3.3/rosetta_source/html/core+numeric+protocols/dc/db0/_atom___8hh-source.html#l00475
      For definition
    */


CL_DEFMETHOD Vector3 Joint_O::position(chem::NVector_sp coords) const
{
  if (this->_PositionIndexX3>=0 && this->_PositionIndexX3 <coords->length()) {
    Vector3 pos(  (*coords)[this->_PositionIndexX3+0],
                  (*coords)[this->_PositionIndexX3+1],
                  (*coords)[this->_PositionIndexX3+2]);
    return pos;
  }
  SIMPLE_ERROR("Out of range Joint_O::position {} at {} when coord length is {} _EndPositionIndexX3 = {}",
               _rep_(this->asSmartPtr()), this->_PositionIndexX3, coords->length(), this->_EndPositionIndexX3);
}

CL_DEFMETHOD void Joint_O::setPosition(chem::NVector_sp coords,const Vector3& pos)
{
  if (this->_PositionIndexX3>=0 && this->_PositionIndexX3 <coords->length()) {
    (*coords)[this->_PositionIndexX3+0] = pos.getX();
    (*coords)[this->_PositionIndexX3+1] = pos.getY();
    (*coords)[this->_PositionIndexX3+2] = pos.getZ();
    return;
  }
  SIMPLE_ERROR("Out of range Joint_O::setPosition {} at {} when coord length is {} _EndPositionIndexX3 = {}",
               _rep_(this->asSmartPtr()), this->_PositionIndexX3, coords->length(), this->_EndPositionIndexX3);
}


CL_LISPIFY_NAME("clearProperty");
CL_DEFMETHOD void Joint_O::clearProperty(core::Symbol_sp symbol)
{
  this->_Properties = core::core__rem_f(this->_Properties,symbol);
}

CL_DOCSTRING(R"dx(Set the property **symbol** of **this** (a chem:matter) to **value**.)dx");
CL_DEFMETHOD void Joint_O::setProperty(core::Symbol_sp symbol, core::T_sp value)
{
  this->_Properties = core::core__put_f(this->_Properties,value,symbol);
}


CL_DOCSTRING(R"dx(Return the property **symbol** of **this** (a chem:matter) - if it isn't defined return NIL.)dx");
CL_LISPIFY_NAME("getProperty");
CL_DEFMETHOD core::T_sp Joint_O::getProperty(core::Symbol_sp symbol)
{
  core::T_sp res = core::cl__getf(this->_Properties,symbol,unbound<core::T_O>());
  if (res.unboundp()) {
    stringstream props;
    props << _rep_(this->_Properties);
    SIMPLE_ERROR("You asked for an unknown property[{}] for matter[{}@{}] - the available properties are[{}]" , _rep_(symbol) , this->__repr__() , (void*)this , props.str()  );
  }
  return res;
}

CL_DOCSTRING(R"dx(Return the property **symbol** of **this** (a chem:matter) - if it isn't defined return **defval**.)dx");
CL_LISPIFY_NAME("getPropertyOrDefault");
CL_DEFMETHOD core::T_sp Joint_O::getPropertyOrDefault(core::Symbol_sp prop,core::T_sp defval)
{
  return core::cl__getf(this->_Properties,prop,defval);
}

CL_DOCSTRING(R"dx(Return T if the property **symbol** of **this** (a chem:matter) is defined.)dx");
CL_LISPIFY_NAME("hasProperty");
CL_DEFMETHOD bool Joint_O::hasProperty(core::Symbol_sp symbol)
{
  return !core::cl__getf(this->_Properties,symbol,unbound<core::T_O>()).unboundp();
}
       
CL_DEFMETHOD void Joint_O::walkChildren(core::Function_sp callback)
{
  LOG("There are {} children" , this->_numberOfChildren() );
  for ( int i=0; i<this->_numberOfChildren(); i++ )
  {
    Joint_sp child = this->_child(i);
    core::eval::funcall(callback,child);
    child->walkChildren(callback);
  }
}

SYMBOL_EXPORT_SC_(KinPkg,joint_calculate_position_index);
int Joint_O::calculatePositionIndex(chem::AtomTable_sp atomTable) {
  core::T_sp index = core::eval::funcall(_sym_joint_calculate_position_index, this->asSmartPtr(), atomTable);
  if (index.fixnump()) {
    int indexX3 = index.unsafe_fixnum()*3;
    return indexX3;
  }
  SIMPLE_ERROR("joint-calculate-position-index did not return a fixnum");
}


CL_DEFUN void kin__walk(Joint_sp joint, core::Function_sp callback) {
  core::eval::funcall(callback, joint);
  joint->walkChildren(callback);
}

void Joint_O::walkResidueTree(int residueId, core::Function_sp callback)
{
  LOG("There are {} children" , this->_numberOfChildren() );
  for ( int i=0; i<this->_numberOfChildren(); i++ )
  {
    Joint_sp child = this->_child(i);
    if ( child->atomId().residueId() != residueId ) continue;
    core::eval::funcall(callback,child);
    child->walkResidueTree(residueId,callback);
  }
}

CL_DEFMETHOD void Joint_O::updateInternalCoord(chem::NVector_sp coords)
{
  this->_updateInternalCoord(coords);
}

CL_DEFMETHOD void Joint_O::updateInternalCoords(chem::NVector_sp coords)
{
  this->_updateInternalCoord(coords);
  for (int childIdx=0; childIdx<this->_numberOfChildren(); childIdx++ ) {
    this->_child(childIdx)->updateInternalCoords(coords);
  }
}

void Joint_O::_updateChildrenXyzCoords(chem::NVector_sp coords) {
  for ( int ii=0; ii < this->_numberOfChildren(); ii++) {
    Stub stub = this->_child(ii)->getInputStub(coords);
    // I should ratchet the newStub around the X axis and use relative dihedral
    this->_child(ii)->_updateXyzCoord(coords,stub);
    // ratchet newStub
//    this->_DofChangePropagatesToYoungerSiblings = false;
    this->noteXyzUpToDate();
  }
  for ( int ii=0; ii < this->_numberOfChildren(); ii++) {
    this->_child(ii)->_updateChildrenXyzCoords(coords);
  }
}


CL_DEFMETHOD void Joint_O::updateXyzCoord(chem::NVector_sp coords)
{
  KIN_LOG("base method\n");
  Stub stub = this->getInputStub(coords);
  this->_updateXyzCoord(coords,stub);
}


CL_DEFMETHOD void Joint_O::updateXyzCoords(chem::NVector_sp coords)
{
  Stub stub = this->getInputStub(coords);
  this->_updateXyzCoord(coords,stub);
  this->_updateChildrenXyzCoords(coords);
}

CL_DEFMETHOD core::T_sp Joint_O::getParentOrNil() const {
  return this->_Parent.unboundp() ? nil<core::T_O>() : gc::As_unsafe<core::T_sp>(this->_Parent);
};


SYMBOL_SC_(KinPkg,atom);
core::Symbol_sp Joint_O::typeSymbol() const { return _sym_atom;};

};
