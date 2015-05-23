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
    EXPOSE_CLASS(kinematics,Atom_O);




    void Atom_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Atom_O>()
	    .def("parent",&Atom_O::parent)
	    .def("typeSymbol",&Atom_O::typeSymbol)
	    .def("asString",&Atom_O::asString)
	    .def("appendChild",&Atom_O::appendChild)
	    .def("insertChild",&Atom_O::insertChild)
//	    .def("setComment",&Atom_O::setComment)
	    .def("atomId",&Atom_O::atomId)
	    .def("getPosition",&Atom_O::getPosition)
	    .def("setPosition",&Atom_O::setPosition)
	    .def("comment",&Atom_O::comment)
	    .def("stubAtom1",&Atom_O::stubAtom1)
	    .def("stubAtom2",&Atom_O::stubAtom2)
	    .def("stubAtom3",&Atom_O::stubAtom3)
	    .def("indexOfChild",&Atom_O::indexOfChild)
	    .def("atom-id-values",&Atom_O::atomIdValues)
	    ;
    }

    void Atom_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,Atom,"","",_lisp)
	    .def("parent",&Atom_O::parent)
	    .def("typeSymbol",&Atom_O::typeSymbol)
	    .def("asString",&Atom_O::asString)
	    .def("appendChild",&Atom_O::appendChild)
	    .def("insertChild",&Atom_O::insertChild)
//	    .def("setComment",&Atom_O::setComment)
	    .def("atomId",&Atom_O::atomId)
	    .def("getPosition",&Atom_O::getPosition)
	    .def("setPosition",&Atom_O::setPosition)
	    .def("comment",&Atom_O::comment)
	    .def("stubAtom1",&Atom_O::stubAtom1)
	    .def("stubAtom2",&Atom_O::stubAtom2)
	    .def("stubAtom3",&Atom_O::stubAtom3)
	    .def("indexOfChild",&Atom_O::indexOfChild)
	    .def("atom-id-values",&Atom_O::atomIdValues)
	    ;
#endif
    }

    Atom_O::Atom_O(const AtomTree_sp& tree, const RefCountedAtomHandle& handle) :
	_AtomTree(tree),
	_Handle(handle)
    {};


	//! Copy ctor increments ref count
    Atom_O::Atom_O(const Atom_O& other) :
	    _AtomTree(other._AtomTree),
	    _Handle(other._Handle)
	    {};

    Atom_O::~Atom_O()
	{
	    if ( this->_AtomTree.pointerp() )
	    {
		this->_Handle.reset();
	    }
	}

    Atom_O& Atom_O::operator=(const Atom_O& other)
	    {_OF();
		ASSERTF(this->_AtomTree == other._AtomTree,BF("You tried to copy Atom_Os that are not referencing the same AtomTree"));
		// Assignment operator
		if ( this != &other ) // avoid self-assignment
		{
		    // If we already point to something then dereference it
		    if ( this->_AtomTree.pointerp() && this->_AtomTree.notnilp() )
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
	    core::Fixnum_O::create(this->get()->id().moleculeId()),
	    core::Fixnum_O::create(this->get()->id().residueId()),
	    core::Fixnum_O::create(this->get()->id().atomId()));
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

    Atom_sp Atom_O::parent() const
    {_OF();
	if ( this->_Handle.notDefined() ) return Atom_O::_nil;
	const Atom* atom = this->_Handle.get();
	return Atom_O::create(atom->parent());
    }

    Atom_sp Atom_O::stubAtom1() const
    {_OF();
	if ( this->_Handle.notDefined() ) return Atom_O::_nil;
	const Atom* atom = this->_Handle.get();
	RefCountedAtomHandle sa = atom->stubAtom1();
	if ( sa.notDefined() ) return Atom_O::_nil;
	return Atom_O::create(sa);
    }

    Atom_sp Atom_O::stubAtom2() const
    {_OF();
	if ( this->_Handle.notDefined() ) return Atom_O::_nil;
	const Atom* atom = this->_Handle.get();
	RefCountedAtomHandle sa = atom->stubAtom2();
	if ( sa.notDefined() ) return Atom_O::_nil;
	return Atom_O::create(sa);
    }

    Atom_sp Atom_O::stubAtom3(AtomTree_sp at) const
    {_OF();
	if ( this->_Handle.notDefined() ) return Atom_O::_nil;
	const Atom* atom = this->_Handle.get();
	RefCountedAtomHandle sa = atom->stubAtom3(at);
	if ( sa.notDefined() ) return Atom_O::_nil;
	return Atom_O::create(sa);
    }


}; /* kinematics */
