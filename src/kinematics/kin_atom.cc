#define DEBUG_LEVEL_FULL


#include <clasp/core/foundation.h>
#include <cando/kinematics/symbolTable.h>
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/atomo.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/atom.h>

namespace kinematics
{

#define	ASSERT_VALID_HANDLE(tree,handle)				\
    ASSERTF((int)handle<tree->numberOfEntries(),BF("The handle[%d] is out of range (0->%d]") % handle % tree->numberOfEntries()); \
    ASSERTF(tree->_AtomHolders[handle]._Type != unused,BF("The handle represents an unused node"));



    core::Lisp_sp Atom::lisp() const
    {
	return this->_Me.lisp();
    }

    string Atom::asString() const
    {_OF();
	stringstream ss;
	ss << "Node[" << this->typeSymbol()->__repr__() << "] RefCount["<<this->_Me.refCount()<<"] id=" << this->_Id.asString() << std::endl;
	ss << "Children: ";
	if ( this->_numberOfChildren() == 0 )
	{
	    ss << "-none-";
	} else
	{
	    for ( int i=0; i<(int)this->_numberOfChildren(); i++ )
	    {
		Atom_sp atom = Atom_O::create(this->_child(i));
		ss << " [" << atom->handle().holderIndex() << "]";
		ss << atom->handle().holder()->typeAsString();
	    }
	}
	ss << std::endl;
	ss << "Position: " << this->_Position.asString() << std::endl;
	return ss.str();
    }


    void Atom::insertChild(int before, RefCountedAtomHandle& child )
    {_OF();
	ASSERTF(child.get() != this,BF("Circular atom reference"));
	this->_insertChild(before,child);
	child.get()->setParent(this->_Me);
    }

    void Atom::appendChild(RefCountedAtomHandle& child)
    {_OF();
	ASSERTF(child.get() != this,BF("Circular atom reference this@%p child@%p")
		% this % (child.get()));
	if ( child.isJump() )
	{
	    int idx = this->firstNonJumpChildIndex();
	    this->_insertChild(idx,child);
	    child.get()->setParent(this->_Me);
	    return;
	}
	this->_appendChild(child);
	child.get()->setParent(this->_Me);
    }



    void Atom::insertChild(RefCountedAtomHandle& child)
    {_OF();
	ASSERTF(child.get() != this,BF("Circular atom reference"));
	LOG(BF("Inserting child of type: %s") % child.holder()->typeAsString() );
	if ( child.isJump() )
	{
	    LOG(BF("It's a jump, inserting it at the start"));
	    this->_insertChild(0,child);
	    child.get()->setParent(this->_Me);
	} else
	{
	    LOG(BF("It's a non-jump atom"));
	    int firstNonJumpIndex = this->firstNonJumpChildIndex();
	    if ( firstNonJumpIndex < this->_numberOfChildren() )
	    {
		IMPLEMENT_ME(); // Really? Is the below just for debugging?
#ifdef DEBUG_ON
		RefCountedAtomHandle& firstNonJumpHandle = this->_child(firstNonJumpIndex);
		LOG(BF("The current firstNonJumpHandle is of type: %s")
		    % firstNonJumpHandle.holder()->typeAsString());
		{
		    LOG(BF("We are inserting a BondedAtom"));
		}
#endif
	    } else
	    {
		LOG(BF("We are at the end of the Children - appending"));
	    }
	    this->insertChild(firstNonJumpIndex,child);
	}
    }

    RootAtomInfo const* Atom::rootAtomInfo() const
    {
	THROW_HARD_ERROR(BF("rootAtomInfo is not available for non-root atoms"));
    }

    int Atom::firstNonJumpChildIndex() const
    {_OF();
	for ( int i=0; i<this->_numberOfChildren(); i++ )
	{
	    if ( !this->_child(i).isJump() ) return i;
	}
	return this->_numberOfChildren();
    }




    void Atom::recursiveDumpChildrenIntoStringStream(const string& prefix, stringstream& out)
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
	    Atom* child = this->_child(i).get();
	    child->recursiveDumpChildrenIntoStringStream(newPrefix,out);
	}
    }


    bool Atom::stubDefined() const
    {
	// Inspired by rosetta::Atom_::stub_defined()
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
	if ( this->isJump() )
	{
	    RefCountedAtomHandle first = this->getNonJumpAtom(0);
	    if ( first.isDefined() &&
		 (first.get()->getNonJumpAtom(0).isDefined() || this->getNonJumpAtom(1).isDefined() ))
	    {
		return true;
	    } else
	    {
		return false;
	    }
	}
	return true;
    }


    RefCountedAtomHandle Atom::getNonJumpAtom(int offset) const
    {
	int idx = this->firstNonJumpChildIndex();
	idx += offset;
	if ( idx >= this->_numberOfChildren() ) return RefCountedAtomHandle::undefined();
	return this->_child(idx);
    }


    RefCountedAtomHandle Atom::previousChild(const RefCountedAtomHandle& ch) const
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
	    return RefCountedAtomHandle::undefined();
	}
	return this->_child(ii-1);
    }


    /*! See file:///Users/meister/Development/rosetta3.3/rosetta_source/html/core+numeric+protocols/dc/db0/_atom___8hh-source.html#l00475
      For definition
    */
    RefCountedAtomHandle Atom::inputStubAtom3(AtomTree_sp at) const
    {_OF();
	ASSERTF(this->parent().isDefined(),BF("The parent isn't defined"));
	RefCountedAtomHandle sibling = this->previousSibling();
	if ( this->isJump() || sibling.notDefined() || sibling.get()->isJump() ||
	     (this->parent().get()->isJump() && sibling.get()->id() == this->parent().get()->stubAtom2Id() ) )
	{
	    return this->parent().get()->stubAtom3(at);
	} else
	{
	    return sibling;
	}
    }



    Stub Atom::getStub(AtomTree_sp at) const
    {_OF();
	LOG(BF("Getting stub for atom[%s]") % this->asString() );
	LOG(BF("stubAtom1(): %s") % this->stubAtom1().get()->asString() );
	LOG(BF("stubAtom2(): %s") % this->stubAtom2().get()->asString() );
	LOG(BF("stubAtom3(): %s") % this->stubAtom3(at).get()->asString() );
	return Stub(
	    this->position(),
	    this->stubAtom1().get()->position(),
	    this->stubAtom2().get()->position(),
	    this->stubAtom3(at).get()->position());
    };


    Stub Atom::getInputStub(AtomTree_sp at) const
    {_OF();
	if ( this->parent().isDefined() )
	{
	    return Stub(
		this->inputStubAtom0().get()->position(),
		this->inputStubAtom1().get()->position(),
		this->inputStubAtom2().get()->position(),
		this->inputStubAtom3(at).get()->position());
	} else
	{
	    return Stub();
	}
    }

    void Atom::walkChildren(const AtomTreeWalkFunctor& functor)
    {_OF();
	LOG(BF("There are %d children") % this->_numberOfChildren() );
	for ( int i=0; i<this->_numberOfChildren(); i++ )
	{
	    RefCountedAtomHandle handle = this->_child(i);
	    Atom_sp child = Atom_O::create(handle);
	    functor(child);
	    child->get()->walkChildren(functor);
	}
    }



    void Atom::walkResidueTree(int residueId, const AtomTreeWalkFunctor& functor)
    {_OF();
	LOG(BF("There are %d children") % this->_numberOfChildren() );
	for ( int i=0; i<this->_numberOfChildren(); i++ )
	{
	    RefCountedAtomHandle handle = this->_child(i);
	    Atom_sp child = Atom_O::create(handle);
	    if ( child->atomId().residueId() != residueId ) continue;
	    functor(child);
	    child->get()->walkResidueTree(residueId,functor);
	}
    }


    void Atom::updateXyzCoords(AtomTree_sp at)
    {_OF();
	Stub stub(this->getInputStub(at));
	this->updateXyzCoords(stub,at);
    }



    SYMBOL_SC_(KinPkg,atom);
    core::Symbol_sp Atom::typeSymbol() const {_OF(); return _sym_atom;};
    



};





