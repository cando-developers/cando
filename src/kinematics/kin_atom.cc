/*
    File: kin_atom.cc
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
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/atom.h>


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
    SIMPLE_ERROR(BF("Cannot convert CoordinateCalculator %p") % v);
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
    core::Symbol_sp sym = o.as<core::Symbol_O>();
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
      SIMPLE_ERROR(BF("Cannot convert CoordinateCalculator %s") % _rep_(o));
    }
  }
};

  

};



namespace kinematics {
FORWARD(JumpAtom);

#define	ASSERT_VALID_HANDLE(tree,handle)				\
  ASSERTF((int)handle<tree->numberOfEntries(),BF("The handle[%d] is out of range (0->%d]") % handle % tree->numberOfEntries()); \
  ASSERTF(tree->_AtomHolders[handle]._Type != unused,BF("The handle represents an unused node"));


CL_DEFMETHOD void Joint_O::setToInternal(core::Symbol_sp cc) {this->_ToInternal = translate::from_object<CoordinateCalculator>(cc)._v; };
CL_DEFMETHOD core::Symbol_sp Joint_O::getToInternal() {return translate::to_object<CoordinateCalculator>::convert(this->_ToInternal); };
CL_DEFMETHOD void Joint_O::setToExternal(core::Symbol_sp cc) {this->_ToExternal = translate::from_object<CoordinateCalculator>(cc)._v; };
CL_DEFMETHOD core::Symbol_sp Joint_O::getToExternal() {return translate::to_object<CoordinateCalculator>::convert(this->_ToExternal); };


string Joint_O::__repr__() const {
  stringstream ss;
  ss << "#<";
  ss << this->className();
  ss << " " << _rep_(this->_Name) << ">";
  return ss.str();
}
 
void Joint_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,parent),this->_Parent);
  node->field(INTERN_(kw,name),this->_Name); // name
  node->field(INTERN_(kw,id),this->_Id);
  node->field(INTERN_(kw,to_external),this->_ToExternal);
  node->field(INTERN_(kw,to_internal),this->_ToInternal);
  node->field_if_not_default(INTERN_(kw,pos),this->_Position, Vector3());
}

CL_DEFMETHOD core::T_sp Joint_O::name() const {
  return this->_Name;
}


string Joint_O::asString() const
{_OF();
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
  ss << "Position: " << this->_Position.asString() << std::endl;
  return ss.str();
}

void Joint_O::setParent(Joint_sp parent)
{
  this->_Parent = parent;
}
    
    

void Joint_O::insertChild(int before, Joint_sp child )
{_OF();
  ASSERTF(child.get() != this,BF("Circular atom reference"));
  this->_insertChild(before,child);
  child.get()->setParent(this->asSmartPtr());
}

void Joint_O::appendChild(Joint_sp child)
{_OF();
  ASSERTF(child.get() != this,BF("Circular atom reference this@%p child@%p")
          % this % (child.get()));
  if ( gc::IsA<JumpJoint_sp>(child) )
  {
    int idx = this->firstNonJumpChildIndex();
    this->_insertChild(idx,child);
    child.get()->setParent(this->asSmartPtr());
    return;
  }
  this->_appendChild(child);
  child.get()->setParent(this->asSmartPtr());
}

void Joint_O::eraseChild(Joint_sp child)
{_OF();
  Joint_sp atom = this->asSmartPtr();
  for (int i=0; i<atom->_numberOfChildren(); i++ ) {
    if (child == atom->_child(i)) {
      atom->_releaseChild(i);
      return;
    }
  }
  SIMPLE_ERROR(BF("Could not find child %s in %s - so it could not be erased") % _rep_(child) % _rep_(this->asSmartPtr()));
}


int Joint_O::indexOfChild(Joint_sp child)
{_OF();
  Joint_sp atom = this->asSmartPtr();
  for (int i=0; i<atom->_numberOfChildren(); i++ )
  {
    if (child == atom->_child(i)) return i;
  }
  SIMPLE_ERROR(BF("Could not find child"));
}




void Joint_O::insertChild(Joint_sp child)
{_OF();
  ASSERTF(child.get() != this,BF("Circular atom reference"));
  LOG(BF("Inserting child: %s") % _rep_(child));
  if ( gc::IsA<JumpJoint_sp>(child) )
  {
    LOG(BF("It's a jump, inserting it at the start"));
    this->_insertChild(0,child);
    child.get()->setParent(this->asSmartPtr());
  } else {
    LOG(BF("It's a non-jump atom"));
    int firstNonJumpIndex = this->firstNonJumpChildIndex();
    if ( firstNonJumpIndex < this->_numberOfChildren() )
    {
      //  IS THIS CORRECT?????
      this->_insertChild(firstNonJumpIndex,child);
#ifdef DEBUG_ON
      Joint_sp firstNonJumpHandle = this->_child(firstNonJumpIndex);
      LOG(BF("The current firstNonJumpHandle is of type: %s")
          % _rep_(firstNonJumpHandle));
      {
        LOG(BF("We are inserting a BondedAtom"));
      }
#endif
    } else
    {
      LOG(BF("We are at the end of the Children - appending"));
      //  IS THIS CORRECT?????
      this->_appendChild(child);
    }
    // this->insertChild(firstNonJumpIndex,child);
  }
}

RootAtomInfo const* Joint_O::rootAtomInfo() const
{
  THROW_HARD_ERROR(BF("rootAtomInfo is not available for non-root atoms"));
}

int Joint_O::firstNonJumpChildIndex() const
{_OF();
  for ( int i=0; i<this->_numberOfChildren(); i++ )
  {
    if ( !gc::IsA<JumpJoint_sp>(this->_child(i)) ) return i;
  }
  return this->_numberOfChildren();
}




void Joint_O::recursiveDumpChildrenIntoStringStream(const string& prefix, stringstream& out)
{_OF();
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
	// note -- in counting dependent atoms, exclude JumpAtoms
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
    Joint_sp first = this->getNonJumpAtom(0);
    if ( first.boundp() &&
         (first->getNonJumpAtom(0) || this->getNonJumpAtom(1).boundp() ))
    {
      return true;
    } else
    {
      return false;
    }
  }
  return true;
}


Joint_sp Joint_O::getNonJumpAtom(int offset) const
{
  int idx = this->firstNonJumpChildIndex();
  idx += offset;
  if ( idx >= this->_numberOfChildren() ) return _Unbound<Joint_O>();
  return this->_child(idx);
}

CL_DEFMETHOD core::List_sp Joint_O::children() const
{
  ql::list l;
  for ( int i(0), iEnd(this->_numberOfChildren()); i<iEnd; ++i ) {
    l << this->_child(i);
  }
  return l.cons();
}
    

Joint_sp Joint_O::previousChild(Joint_sp ch) const
{_OF();
  int num = this->_numberOfChildren();
  int ii;
  for ( ii=0; ii<num; ii++ )
  {
    if ( this->_child(ii) == ch ) break;
  }
  if ( ii == num )
  {
    SIMPLE_ERROR(BF("Could not find child"));
  }
  if ( ii == 0 )
  {
    return _Unbound<Joint_O>();
  }
  return this->_child(ii-1);
}


    /*! See file:///Users/meister/Development/rosetta3.3/rosetta_source/html/core+numeric+protocols/dc/db0/_atom___8hh-source.html#l00475
      For definition
    */
Joint_sp Joint_O::inputStubAtom3(AtomTree_sp at) const
{_OF();
  ASSERTF(this->parent().boundp(),BF("The parent isn't defined"));
  if (this->parent().unboundp()) {
    SIMPLE_ERROR(BF("inputStubAtom2 parent of %s isn't defined") % _rep_(this->asSmartPtr()));
  }
  Joint_sp sibling = this->previousSibling();
  if ( gc::IsA<JumpJoint_sp>(this->asSmartPtr())
       || sibling.unboundp()
       || gc::IsA<JumpJoint_sp>(sibling)
       || (gc::IsA<JumpJoint_sp>(this->parent())
           && sibling->id() == this->parent()->stubAtom2Id() ) ) {
    return this->parent()->stubAtom3(at);
  } else {
    return sibling;
  }
}

CL_DEFMETHOD Vector3 Joint_O::getPosition() const
{
  Vector3 pos = this->_Position;
  return pos;
}

CL_DEFMETHOD void Joint_O::setPosition(const Vector3& pos)
{
  this->_Position = pos;
}


CL_LISPIFY_NAME("clearProperty");
CL_DEFMETHOD void Joint_O::clearProperty(core::Symbol_sp symbol)
{
  this->_Properties = core::core__rem_f(this->_Properties,symbol);
}

CL_DOCSTRING("Set the property **symbol** of **this** (a chem:matter) to **value**.");
CL_DEFMETHOD void Joint_O::setProperty(core::Symbol_sp symbol, core::T_sp value)
{
  this->_Properties = core::core__put_f(this->_Properties,value,symbol);
}


CL_DOCSTRING("Return the property **symbol** of **this** (a chem:matter) - if it isn't defined return NIL.");
CL_LISPIFY_NAME("getProperty");
CL_DEFMETHOD core::T_sp Joint_O::getProperty(core::Symbol_sp symbol)
{
  core::T_sp res = core::cl__getf(this->_Properties,symbol,_Unbound<core::T_O>());
  if (res.unboundp()) {
    stringstream props;
    props << _rep_(this->_Properties);
    SIMPLE_ERROR(BF("You asked for an unknown property[%s] for matter[%s@%p] - the available properties are[%s]") % _rep_(symbol) % this->__repr__() % this % props.str()  );
  }
  return res;
}

CL_DOCSTRING("Return the property **symbol** of **this** (a chem:matter) - if it isn't defined return **defval**.");
CL_LISPIFY_NAME("getPropertyOrDefault");
CL_DEFMETHOD core::T_sp Joint_O::getPropertyOrDefault(core::Symbol_sp prop,core::T_sp defval)
{
  return core::cl__getf(this->_Properties,prop,defval);
}

CL_DOCSTRING("Return T if the property **symbol** of **this** (a chem:matter) is defined.");
CL_LISPIFY_NAME("hasProperty");
CL_DEFMETHOD bool Joint_O::hasProperty(core::Symbol_sp symbol)
{
  return !core::cl__getf(this->_Properties,symbol,_Unbound<core::T_O>()).unboundp();
}



       
void Joint_O::walkChildren(core::Function_sp callback)
{_OF();
  LOG(BF("There are %d children") % this->_numberOfChildren() );
  for ( int i=0; i<this->_numberOfChildren(); i++ )
  {
    Joint_sp child = this->_child(i);
    core::eval::funcall(callback,child);
    child->walkChildren(callback);
  }
}



void Joint_O::walkResidueTree(int residueId, core::Function_sp callback)
{_OF();
  LOG(BF("There are %d children") % this->_numberOfChildren() );
  for ( int i=0; i<this->_numberOfChildren(); i++ )
  {
    Joint_sp child = this->_child(i);
    if ( child->atomId().residueId() != residueId ) continue;
    core::eval::funcall(callback,child);
    child->walkResidueTree(residueId,callback);
  }
}

CL_DEFMETHOD void Joint_O::updateInternalCoord()
{
  this->_updateInternalCoord();
}

CL_DEFMETHOD void Joint_O::updateXyzCoord()
{_OF();
  KIN_LOG(BF("base method\n"));
  Stub stub = this->parent()->getStub();
  this->_updateXyzCoord(stub);
}


CL_DEFMETHOD void Joint_O::updateXyzCoords()
{_OF();
  Stub stub = this->parent()->getStub();
  this->_updateXyzCoords(stub);
}

SYMBOL_SC_(KinPkg,atom);
core::Symbol_sp Joint_O::typeSymbol() const {_OF(); return _sym_atom;};

CL_DEFMETHOD core::T_sp Joint_O::getParent() const {
  return this->_Parent.unboundp() ? _Nil<core::T_O>() : gc::As_unsafe<core::T_sp>(this->_Parent);
};


};
