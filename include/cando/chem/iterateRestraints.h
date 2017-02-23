/*
    File: iterateRestraints.h
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

#ifndef	IterateRestraints_H //[
#define IterateRestraints_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/iterator.h>
#include <cando/chem/iterateMatter.h>



#include <cando/chem/chemPackage.h>


namespace chem {


  SMART(IterateRestraints );
  class IterateRestraints_O : public core::Iterator_O
  {
    LISP_CLASS(chem,ChemPkg,IterateRestraints_O,"IterateRestraints",core::Iterator_O);
  public:
    void	initialize();
  private:
    Matter_sp			_Top;
    core::List_sp				_Restraints;
    core::List_sp				_CurRestraint;
  private:
    void _accumulateRestraints(Matter_sp m);

  public:
    static IterateRestraints_sp create(Matter_sp matter);
  public:
    void initTop(Matter_sp matter);
    void first();
    void next();
    bool isDone();
    core::T_sp currentObject();
	
//	IterateRestraints_O( const IterateRestraints_O& ss ); //!< Copy constructor

  IterateRestraints_O() : _Top(_Unbound<Matter_O>()), _Restraints(_Nil<T_O>()), _CurRestraint(_Nil<T_O>()) {};
  };



};
TRANSLATE(chem::IterateRestraints_O);
#endif //]
