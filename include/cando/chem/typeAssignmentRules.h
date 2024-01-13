/*
    File: typeAssignmentRules.h
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
#if 0       
       
#ifndef	TypeAssignmentRules_H //[
#define TypeAssignmentRules_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/chemInfo.h>

#include <cando/chem/chemPackage.h>


namespace chem {

  SMART(OneTypeRule );
  class OneTypeRule_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,OneTypeRule_O,"OneTypeRule",core::CxxObject_O);
  public:
    void	initialize();
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    core::Symbol_sp _AssignTypeName;
    ChemInfo_sp	_Match;
  public:
    void setAssignTypeName(core::Symbol_sp nm) { this->_AssignTypeName = nm;};
    core::Symbol_sp  getAssignTypeName()	{return this->_AssignTypeName;};
    void setCode(ChemInfo_sp code) { this->_Match = code; };

    bool testMatch(Atom_sp a);


    OneTypeRule_O( const OneTypeRule_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(OneTypeRule_O);
  };

  SMART(TypeAssignmentRules );
  class TypeAssignmentRules_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,TypeAssignmentRules_O,"TypeAssignmentRules",core::CxxObject_O);
  public:
    void	initialize();
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  private:
    WildElementDict_sp	_WildElementDict;
    core::MDArrayT_sp _Rules;
  public:
  public:

    CL_DEFMETHOD void setWildElementDict(WildElementDict_sp dict) { this->_WildElementDict = dict;};
    WildElementDict_sp getWildElementDict(WildElementDict_sp dict) { return this->_WildElementDict;};

    CL_DEFMETHOD void appendRule(OneTypeRule_sp rule) { this->_Rules->vectorPushExtend(rule);};

    core::Symbol_sp calculateType(Atom_sp a);

    TypeAssignmentRules_O( const TypeAssignmentRules_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(TypeAssignmentRules_O);
  };



};
#endif //]
#endif
