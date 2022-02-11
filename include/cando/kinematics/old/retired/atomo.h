/*
    File: atomo.h
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
#ifndef	_kinematics_AtomO_H
#define _kinematics_AtomO_H

#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <clasp/core/multipleValues.fwd.h>
#include <cando/kinematics/point.fwd.h>
#include <cando/kinematics/chainNode.fwd.h>
#include <cando/kinematics/pointTree.fwd.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{

    FORWARD(PointTree);

    /*! @class Wraps an PointTree and an AtomHandle and provides
      safe, reference counted access to the underlying Atom.
    */
    class Point_O : public core::CxxObject_O
    {
	LISP_CLASS(kinematics,KinPkg,Point_O,"Atom",core::CxxObject_O);
    private:
	PointTree_sp		_PointTree;
	RefCountedAtomHandle	_Handle;
    public:
    Point_O() : _Handle() {}
	Point_O(const PointTree_sp& tree, Point_sp handle);

	//! Copy ctor increments ref count
	Point_O(const Point_O& other);
	virtual ~Point_O();

    public:
	static Point_sp create(Point_sp handle)
	{
          auto atom = gctools::GC<Point_O>::allocate_with_default_constructor();
	    atom->_PointTree = handle.pointTree();
	    atom->_Handle = handle;
	    return atom;
	}
	static Point_sp create(const PointTree_sp& tree, uint handleIndex)
	{
          auto atom = gctools::GC<Point_O>::allocate_with_default_constructor();
	    atom->_PointTree = tree;
	    RefCountedAtomHandle handle(tree.get(),handleIndex);
	    try
	    {
		atom->_Handle = handle;
		handle.reset();
	    } catch (...)
	      {
		  handle.reset();
		  throw;
	      }
	    return atom;
	}


	RefCountedAtomHandle& handle() { return this->_Handle;};
	RefCountedAtomHandle const & handle() const { return this->_Handle;};


	Atom const& operator*() const
	{
	    return *(this->_Handle.get());
	}

	Atom* operator->()
	{
	    return this->_Handle.get();
	}


	Atom* get()
	{
	    if ( this->_Handle.notDefined()) return NULL;
	    return this->_Handle.get();
	}

	Atom const* get() const
	{
	    if ( this->_Handle.notDefined()) return NULL;
	    return this->_Handle.get();
	}


	Point_O& operator=(const Point_O& other);

	string asString() const;

	/*! Return the type of atom as a string */
	core::Symbol_sp typeSymbol() const;

	/*! Return the AtomId */
	chem::AtomId atomId() const;

	/*! Return the AtomId as MultipleValues */
	core::T_mv atomIdValues() const;

	/*! Return the comment for the Atom */
	string comment() const;

	/*! Append the atom as a child of this atom, this wraps appendChild */
	void appendChild(Point_sp atom);


	/*! Insert the atom as a child of this atom, this wraps insertChild */
	void insertChild(Point_sp atom);


	/*! Erase the atom */
	void eraseChild(Point_sp atom);

	/*! Wrap setComment */
	Point_sp setComment(const string& comment);

	/*! Return the index of the child of this Atom */
	int indexOfChild(Point_sp atom);

	/*! Return a string representation */
	string __repr__() const;

	/*! Set the position of the Atom */
	void setPosition(const Vector3& pos);

	/*! Get the position of the Atom */
	Vector3 getPosition() const;

	/*! Get parent */
        gc::Nilable<Point_sp> parent() const;


	/*! Get stubPoint1 */
        gc::Nilable<Point_sp> stubPoint1() const;

	/*! Get stubPoint2 */
        gc::Nilable<Point_sp> stubPoint2() const;

	/*! Get stubPoint3 */
        gc::Nilable<Point_sp> stubPoint3(PointTree_sp tree) const;


#if 0
	/*! Build a Molecule using the graph of monomers described by MonomerNode as a plan */
	void buildChainUsingchainNode(int moleculeId, ChainNode_sp sequenceNode );
#endif
    };







}; /* kinematics */

TRANSLATE(kinematics::Point_O);

#endif /* _kinematics_Point_H */
