/*
    File: candoDatabaseReference.h
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
#ifndef	CandoDatabaseReference_H //[
#define CandoDatabaseReference_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/entity.h>

#include <cando/chem/chemPackage.h>

//#include "core/objRef.fwd.h"// candoDatabaseReference.h wants ObjRef needs objRef.fwd.h

namespace chem {

SMART(ObjRef);

SMART(CandoDatabaseReference);
class CandoDatabaseReference_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,CandoDatabaseReference_O,"CandoDatabaseReference",core::CxxObject_O);
public:
	void	initialize();
private:
	core::ObjRef_sp	_Reference;
public:
    static CandoDatabaseReference_sp create(core::Lisp_sp e,CandoDatabase_sp bdb, const string& link);
public:


		/*!Lookup and return the object in the database referenced by this one
		 * and cast it to the desired type
		 */
	template <class otype>
          gctools::smart_ptr<otype> get() { return this->oget().as<otype>();};

		//!Lookup and return the object in the database referenced by this one
	core::T_sp oget();


	CandoDatabaseReference_O( const CandoDatabaseReference_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CandoDatabaseReference_O);
};


// Constructor function

};
TRANSLATE(chem::CandoDatabaseReference_O);
#endif //]
