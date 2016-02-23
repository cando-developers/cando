/*
    File: atomo.cc
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
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/geom/vector3.h>
#include <cando/geom/ovector3.h>
#include <cando/kinematics/atomo.h>
#include <clasp/core/multipleValues.h>
#include <cando/kinematics/atomTree.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{






    Atom_O::Atom_O(const AtomTree_sp& tree, const RefCountedAtomHandle& handle) :
	_AtomTree(tree),
	_Handle(handle)
    {};


	//! Copy ctor increments ref count
    Atom_O::Atom_O(const Atom_O& other) :
	    _AtomTree(other._AtomTree),
	    _Handle(other._Handle)
	    {};

    Atom_O::~Atom_O() {}

    Atom_O& Atom_O::operator=(const Atom_O& other)
	    {_OF();
		ASSERTF(this->_AtomTree == other._AtomTree,BF("You tried to copy Atom_Os that are not referencing the same AtomTree"));
		// Assignment operator
		if ( this != &other ) // avoid self-assignment
		{
		    // If we already point to something then dereference it
		    if ( this->_AtomTree.generalp() && this->_AtomTree.notnilp() )
		    {
			this->_Handle = other._Handle;
		    }
		}
		return *this;
	    }


    chem::AtomId Atom_O::atomId() const
    {
	return this->get()->id();
    };

    core::T_mv Atom_O::atomIdValues() const
    {_OF();
	return Values(
	    core::clasp_make_fixnum(this->get()->id().moleculeId()),
	    core::clasp_make_fixnum(this->get()->id().residueId()),
	    core::clasp_make_fixnum(this->get()->id().atomId()));
    }

    core::Symbol_sp Atom_O::typeSymbol() const
    {_OF();
	Atom const* atom = this->get();
	if ( atom == NULL ) return _Nil<core::Symbol_O>();
	return atom->typeSymbol();
    }

    string Atom_O::asString() const
    {
	Atom const* atom = this->get();
	if ( atom == NULL ) return "Undefined";
	return atom->asString();
    }


    string Atom_O::comment() const
    {
	if ( this->_Handle.notDefined() ) return "";
        const Atom* atom = this->_Handle.get();
	return atom->_Comment;
    }

    void Atom_O::appendChild(Atom_sp child)
    {_OF();
	ASSERT(this->_AtomTree == child->_AtomTree);
	Atom* parent = this->get();
	parent->appendChild(child->handle());
    }
	

    void Atom_O::insertChild(Atom_sp child)
    {_OF();
	ASSERT(this->_AtomTree == child->_AtomTree);
	Atom* parent = this->get();
	parent->insertChild(child->handle());
    }
	

    void Atom_O::eraseChild(Atom_sp child)
    {_OF();
	Atom* atom = this->handle().get();
	for (int i=0; i<atom->_numberOfChildren(); i++ )
	{
	    if (child->_Handle.holderIndex() == atom->_child(i).holderIndex())
	    {
		atom->_releaseChild(i);
		return;
	    }
	}
	SIMPLE_ERROR(BF("Could not find child"));
    }


    int Atom_O::indexOfChild(Atom_sp child)
    {_OF();
	Atom* atom = this->handle().get();
	for (int i=0; i<atom->_numberOfChildren(); i++ )
	{
	    if (child->_Handle.holderIndex() == atom->_child(i).holderIndex())
	    {
		return i;
	    }
	}
	SIMPLE_ERROR(BF("Could not find child"));
    }


    string Atom_O::__repr__() const
    {
	stringstream ss;
	if ( this->get() == NULL )
	{
	    return "Atom_O::undefined";
	}
	ss << (BF("%s@%p handle[%d]")
	       % this->get()->typeSymbol()
	       % (void*)(this->get())
	       % this->handle().holderIndex() ).str();
	return ss.str();
    }

    void Atom_O::setPosition(const Vector3& pos)
    {
	if ( this->_Handle.notDefined() ) return;
	this->get()->position(pos);
    }

    Vector3 Atom_O::getPosition() const
    {
	return this->get()->position();
    }

gc::Nilable<Atom_sp> Atom_O::parent() const
    {_OF();
      if ( this->_Handle.notDefined() ) return _Nil<core::T_O>();
	const Atom* atom = this->_Handle.get();
	return Atom_O::create(atom->parent());
    }

gc::Nilable<Atom_sp> Atom_O::stubAtom1() const
    {_OF();
      if ( this->_Handle.notDefined() ) return _Nil<core::T_O>();
	const Atom* atom = this->_Handle.get();
	RefCountedAtomHandle sa = atom->stubAtom1();
	if ( sa.notDefined() ) return _Nil<core::T_O>();
	return Atom_O::create(sa);
    }

gc::Nilable<Atom_sp> Atom_O::stubAtom2() const
{_OF();
  if ( this->_Handle.notDefined() ) return _Nil<core::T_O>();
  const Atom* atom = this->_Handle.get();
  RefCountedAtomHandle sa = atom->stubAtom2();
  if ( sa.notDefined() ) return _Nil<core::T_O>();
  return Atom_O::create(sa);
}

gc::Nilable<Atom_sp> Atom_O::stubAtom3(AtomTree_sp at) const
{_OF();
  if ( this->_Handle.notDefined() ) return _Nil<core::T_O>();
  const Atom* atom = this->_Handle.get();
  RefCountedAtomHandle sa = atom->stubAtom3(at);
  if ( sa.notDefined() ) return _Nil<core::T_O>();
  return Atom_O::create(sa);
}


}; /* kinematics */
