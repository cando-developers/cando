/*
    File: energyComponent.fwd.h
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
#ifndef energyComponent_fwd_H
#define energyComponent_fwd_H


typedef float REAL;
typedef int   INT;

struct SetupAccumulator {

  core::List_sp _Plist;
  core::List_sp _ValidSymbols = nil<core::T_O>();

  SetupAccumulator(core::T_sp className, core::T_sp setup );
  ~SetupAccumulator();


  template <typename SlotType>
  void maybe_apply(core::Symbol_sp name, SlotType& value) __attribute__((optnone)) {
    core::T_sp tvalue = cl__getf(this->_Plist,name,nil<core::T_O>());
    printf("%s:%d:%s %s\n", __FILE__, __LINE__, __FUNCTION__, fmt::format("Searched for {} in {} found: {}\n", _rep_(name), _rep_(this->_Plist), _rep_(tvalue)).c_str());
    if (gc::IsA<core::Number_sp>(tvalue)) {
      double dval = translate::make_from_object<SlotType>(tvalue);
      printf("%s\n", fmt::format("Assigning parameter {} with value {}\n", _rep_(name), dval ).c_str());
      value = dval;
    }
    this->_ValidSymbols = core::Cons_O::create(name,this->_ValidSymbols);
  }
};


namespace  chem
{
FORWARD(EnergyComponent);


}
#endif
