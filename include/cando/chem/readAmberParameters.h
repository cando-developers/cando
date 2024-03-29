
/*
    File: readAmberParameters.h
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
#ifndef	ReadAmberParameters_H //[
#define ReadAmberParameters_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>


namespace chem {

SMART(ForceField);
SMART(FFTypesDb);
SMART(FFNonbondDb);
SMART(FFStretchDb);
SMART(FFAngleDb);
SMART(FFPtorDb);
SMART(FFItorDb);

SMART(ReadAmberParameters );
class ReadAmberParameters_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,ReadAmberParameters_O,"ReadAmberParameters",core::CxxObject_O);
 private:
  FFTypesDb_sp	_Types;
  ForceField_sp	_ForceField;
 private:
  core::T_mv determineParmSetFrcModType(core::T_sp fin);
  FFTypesDb_sp parseTypeRules(core::T_sp fin);
  FFNonbondDb_sp parseMasses(core::T_sp fin, FFNonbondDb_sp);
  FFStretchDb_sp parseStretchDb(core::T_sp fin);
  FFAngleDb_sp parseAngleDb(core::T_sp fin);
  FFPtorDb_sp parsePtorDb(core::T_sp fin, core::T_sp system);
  FFItorDb_sp parseItorDb(core::T_sp fin);
  string parseNonbondLabelKindNB(core::T_sp fin);
  void parseNonbondDb(core::T_sp fin, FFNonbondDb_sp );
  ForceField_sp parseAmberFormattedForceField(core::T_sp strm, core::T_sp system);
  ForceField_sp parseFrcModFile(core::T_sp strm, core::T_sp system);
  ForceField_sp parseAmberFormattedForceFieldOrFrcMod(core::T_sp strm, core::T_sp system);
  void parseAtomEquivalences(core::T_sp fin, FFNonbondDb_sp);
  
 public:

  void readTypes(core::T_sp fin);
  void readParameters(core::T_sp fileName, core::T_sp system);

  CL_LISPIFY_NAME("getTypeRules");
  CL_DEFMETHOD     FFTypesDb_sp getTypeRules() { return this->_Types; };
  CL_DEFMETHOD     void setTypeRules(FFTypesDb_sp types) { this->_Types = types; };
  ForceField_sp getForceField();

  ReadAmberParameters_O( const ReadAmberParameters_O& ss ); //!< Copy constructor

 ReadAmberParameters_O() : _Types(unbound<FFTypesDb_O>()), _ForceField(unbound<ForceField_O>()) {};
};


};
#endif //]
