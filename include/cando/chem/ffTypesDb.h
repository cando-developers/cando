/*
    File: ffTypesDb.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffTypes.h
 *
 *	Maintains databases and structures to store types and type assignement
 *	rules.
 */

#ifndef FFTYPES_H
#define	FFTYPES_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/chemInfo.h>
#include <cando/chem/ffBaseDb.h>


#include <cando/chem/chemPackage.h>


namespace chem {

inline	string	XmlTag_FFTypesDb()	{ return "FFTypesDb";};


///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
class	FFTypesDb_O;
    typedef	gctools::smart_ptr<FFTypesDb_O>	FFTypesDb_sp;
SMART(FFTypesDb);
class FFTypesDb_O : public FFBaseDb_O
{
    LISP_CLASS(chem,ChemPkg,FFTypesDb_O,"FFTypesDb",FFBaseDb_O);

public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
public:
    gctools::Vec0<ChemInfo_sp>	_TypeAssignmentRules;
    /*! New rules take precedence over old rules and the first rule that matches
        is the one that is used, so new rules need to be prepended to the array. */
    void forceFieldMerge(FFBaseDb_sp other) {
      FFTypesDb_sp other_types_db = gc::As<FFTypesDb_sp>(other);
      size_t new_size = other_types_db->_TypeAssignmentRules.size() + this->_TypeAssignmentRules.size();
      gctools::Vec0<ChemInfo_sp> newRules;
      newRules.reserve(new_size);
      size_t idx;
      for ( auto it : other_types_db->_TypeAssignmentRules ) {
        newRules.push_back(it);
      }
      for ( auto it : this->_TypeAssignmentRules ) {
        newRules.push_back(it);
      }
      this->_TypeAssignmentRules.swap(newRules);
    }

    CL_LISPIFY_NAME("FFTypesDb-add");
    CL_DEFMETHOD void	add( chem::ChemInfo_sp ci ) {
	this->_TypeAssignmentRules.push_back(ci);
    }

CL_LISPIFY_NAME("FFTypesDb-numberOfRules");
CL_DEFMETHOD     int	numberOfRules() { return this->_TypeAssignmentRules.size();};
    chem::ChemInfo_sp	getRule(uint index);

    void    assignTypes( chem::Matter_sp matter );
    core::Symbol_sp    assignType( chem::Atom_sp atom, bool verbose=false );
    void	initialize();

    DEFAULT_CTOR_DTOR(FFTypesDb_O);
};




};

#endif
