/*
    File: ffBaseDb.h
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
 *	ffBase.h
 *
 *	Maintains a database of stretch types
 */

#ifndef FFBASE_H
#define	FFBASE_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>


#include <cando/chem/chemPackage.h>


namespace chem {

SMART(ForceField);

};

namespace chem {
  typedef	enum { parameterized, estimated, rough, unknown } ParameterizationLevel;

};


namespace chemkw {
    extern core::Symbol_sp& _sym_parameterized;
    extern core::Symbol_sp& _sym_estimated;
    extern core::Symbol_sp& _sym_rough;
    extern core::Symbol_sp& _sym_unknown;
};

namespace translate {
    template <> struct to_object<chem::ParameterizationLevel>
    {								
	typedef	chem::ParameterizationLevel	GivenType;	
	static core::T_sp convert(const GivenType& val)
	{_G();
	    SYMBOL_EXPORT_SC_(ChemKwPkg,parameterized);
	    SYMBOL_EXPORT_SC_(ChemKwPkg,estimated);
	    SYMBOL_EXPORT_SC_(ChemKwPkg,rough);
	    switch (val) {
	    case chem::parameterized:
		return (chemkw::_sym_parameterized);
	    case chem::estimated:
		return (chemkw::_sym_estimated);
	    case chem::rough:
		return (chemkw::_sym_rough);
            case chem::unknown:
                return (chemkw::_sym_unknown);
	    };
	    SIMPLE_ERROR("Illegal ParameterizationLevel");
	}
    };

    template <>
    struct from_object<chem::ParameterizationLevel>
    {								
	typedef	chem::ParameterizationLevel 	ExpectedType;
	typedef	chem::ParameterizationLevel 	DeclareType;
	DeclareType _v;
	from_object(gctools::smart_ptr<core::T_O> o) {
	    if ( o == chemkw::_sym_parameterized ) {
		_v = chem::parameterized;
	    } else if ( o == chemkw::_sym_estimated ) {
		_v = chem::estimated;
	    } else if ( o == chemkw::_sym_rough ) {
		_v = chem::rough;
	    } else {
		SIMPLE_ERROR("Illegal value for ParameterizationLevel {}" , _rep_(o) );
	    }
	}
    };


};


namespace chem {




typedef	enum { stretch, angle, ptor, itor, nonbond, vdw } ParameterType;


class FFParameter_O;
    typedef	gctools::smart_ptr<FFParameter_O>	FFParameter_sp;
SMART(FFParameter);
class FFParameter_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,FFParameter_O,"FFParameter",core::CxxObject_O);

public:
  FFParameter_O() : _Comment(nil<core::T_O>()) {};
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
  ParameterizationLevel	_Level;
  core::T_sp _Comment;

  core::T_sp getComment() const;
  void setComment(core::T_sp comment);

/*! Return the parameterization level of this parameter.
 * In most cases if it isn't "excellent" then it will be noted and the
 * user can go and calculate a good parameter and put it into the force field.
 */
  ParameterizationLevel	level() { return this->_Level;};

/*! Describe this parameter in terms of its ParameterizationLevel and the atom types involved.
 * This will be used when collecting missing or poorly parameterized parameters.
 */
virtual	string levelDescription();
  virtual	ParameterType type() {SUBIMP();};

};



SMART(FFBaseDb);
class FFBaseDb_O : public core::CxxObject_O {
  LISP_CLASS(chem,ChemPkg,FFBaseDb_O,"FFBaseDb",core::CxxObject_O);
 public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
  void initialize();
 public:
  virtual void forceFieldMerge(FFBaseDb_sp other);
  FFBaseDb_O() {};
};

 SMART(FFParameterBaseDb);
class FFParameterBaseDb_O : public FFBaseDb_O {
  LISP_CLASS(chem,ChemPkg,FFParameterBaseDb_O,"FFParameterBaseDb",FFBaseDb_O);
 public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
  void initialize();
 public:
  core::HashTable_sp _Parameters;
  core::ComplexVector_T_sp _ParameterVector;
 public:
  virtual void forceFieldMerge(FFBaseDb_sp other);
 FFParameterBaseDb_O() : FFBaseDb_O(), _Parameters(unbound<core::HashTable_O>()) {};
  core::HashTable_sp parameters() const;
  string __repr__() const;
};


};
#endif
