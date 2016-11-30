/*
    File: chiDihedrals.cc
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
#include <clasp/core/environment.h>
#include <cando/kinematics/chiDihedrals.h>
#include <clasp/core/multipleValues.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//





#if INIT_TO_FACTORIES

#define ARGS_ChiDihedral_O_make "(atom1-index atom2-index atom3-index atom4-index)"
#define DECL_ChiDihedral_O_make ""
#define DOCS_ChiDihedral_O_make "make ChiDihedral"
  ChiDihedral_sp ChiDihedral_O::make(const int atom1Index, const int atom2Index, const int atom3Index, const int atom4Index)
    {
	GC_ALLOCATE(ChiDihedral_O, me );
      me->_Atom1 = atom1Index;
      me->_Atom2 = atom2Index;
      me->_Atom3 = atom3Index;
      me->_Atom4 = atom4Index;
      return me;
    };

#else
    core::T_sp ChiDihedral_O::__init__(::core::Function_sp exec, ::core::Cons_sp args, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {
	this->_Atom1 = env->lookup(Pkg(),"atom1Index").as<core::Rational_O>()->as_int();
	this->_Atom2 = env->lookup(Pkg(),"atom2Index").as<core::Rational_O>()->as_int();
	this->_Atom3 = env->lookup(Pkg(),"atom3Index").as<core::Rational_O>()->as_int();
	this->_Atom4 = env->lookup(Pkg(),"atom4Index").as<core::Rational_O>()->as_int();
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void ChiDihedral_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	node->attribute("a1",this->_Atom1);
	node->attribute("a2",this->_Atom2);
	node->attribute("a3",this->_Atom3);
	node->attribute("a4",this->_Atom4);
    }
#endif


    void ChiDihedral_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    core::T_mv ChiDihedral_O::atomIndices() const
    {_OF();
	return Values(
	    core::clasp_make_fixnum(this->_Atom1),
	    core::clasp_make_fixnum(this->_Atom2),
	    core::clasp_make_fixnum(this->_Atom3),
	    core::clasp_make_fixnum(this->_Atom4));
    }
    







#if INIT_TO_FACTORIES
#else
    core::T_sp ChiList_O::__init__(::core::Function_sp exec, ::core::Cons_sp args, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void ChiList_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->archiveVector0("chis",this->_Chis);
	// Archive other instance variables here
    }
#endif


    void ChiList_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    void ChiList_O::append(ChiDihedral_sp chi)
    {_OF();
	this->_Chis.push_back(chi);
    }

    

}; /* kinematics */
