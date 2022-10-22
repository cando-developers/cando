/*
    File: angle.cc
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
#include <cando/chem/atom.h>
#include <cando/chem/angle.h>
#include <clasp/core/wrappers.h>
namespace chem
{

// ----------------------------------------------------------------------
//








    Angle_sp Angle_O::create(Atom_sp a1, Atom_sp a2, Atom_sp a3)
    {
      auto  obj  = gctools::GC<Angle_O>::allocate_with_default_constructor();
	obj->setup(a1,a2,a3);
	return obj;
    }


    void Angle_O::setup(Atom_sp a1, Atom_sp a2, Atom_sp a3)
    {
	this->_a1 = a1;
	this->_a2 = a2;
	this->_a3 = a3;
    }


#ifdef XML_ARCHIVE
    void Angle_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void Angle_O::initialize()
    {
        this->Base::initialize();
    }


    string Angle_O::__repr__() const
    {
	stringstream ss;
	ss << this->className() << "["<<this->_a1->getName()<<"-"<<this->_a2->getName()<<"-"<<this->_a3->getName()<<"]";
	return ss.str();
    }



}; /* chem */
