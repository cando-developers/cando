/*
    File: bondId.h
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
#ifndef	_BONDID_H
#define _BONDID_H

#include <clasp/core/common.h>
#include <cando/kinematics/jointTree.fwd.h>
#include <cando/kinematics/bondId.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>



namespace kinematics
{
  class BondId_O : public core::CxxObject_O
  {
    LISP_CLASS(kinematics,KinPkg,BondId_O,"BondId",core::CxxObject_O);
  public:
  BondId_O() : _Parent(unbound<Joint_O>()), _Child(unbound<Joint_O>()) {};
  public: // instance variables here
    Joint_sp	_Parent;
    Joint_sp	_Child;
  public:
    static BondId_sp create(Joint_sp parent, Joint_sp child = unbound<Joint_O>());
    string __repr__() const;
  };

}; /* kinematics */

#endif /* _BONDID_H */
