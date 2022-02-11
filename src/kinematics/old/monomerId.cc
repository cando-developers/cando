/*
    File: monomerId.cc
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
#define DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <clasp/core/lispStream.h>
#include <cando/kinematics/jointTree.h>
#include <cando/kinematics/monomerId.h>
#include <clasp/core/wrappers.h>



namespace kinematics
{

string MonomerId::asString() const
{
  stringstream ss;
  ss << "Chain["<<this->_Chain;
  ss << "]Res["<<this->_Monomer <<"]";
  return ss.str();
}


DOCGROUP(cando)
CL_DEFUN void kin__dumpMonomerId(const kinematics::MonomerId& id) {
  core::write_bf_stream(BF("MonomerId %lu %lu\n") % id._Chain % id._Monomer );
}





#if INIT_TO_FACTORIES

CL_LISPIFY_NAME(make-MonomerId);
MonomerId_sp MonomerId_O::make(const int chain, const int monomer)
{
  auto  me  = gctools::GC<MonomerId_O>::allocate_with_default_constructor();
  me->_MonomerId._Chain = chain;
  me->_MonomerId._Monomer = monomer;
  return me;
};

#else
core::T_sp MonomerId_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp)
{
  this->_MonomerId._Chain = translate::from_object<int>::convert(env->lookup(Pkg(),"chain"));
  this->_MonomerId._Monomer = translate::from_object<int>::convert(env->lookup(Pkg(),"monomer"));
  return _Nil<core::T_O>();
};
#endif


#ifdef XML_ARCHIVE
void MonomerId_O::archiveBase(core::ArchiveP node)
{
  node->attribute("Chain",this->_MonomerId._Chain);
  node->attribute("Monomer",this->_MonomerId._Monomer);
}
#endif

};
