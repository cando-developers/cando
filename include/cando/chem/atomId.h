/*
    File: atomId.h
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
#ifndef	chem_atomId_H
#define chem_atomId_H

#include <clasp/core/common.h>
#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/chem/chemPackage.h>


namespace chem
{

class AtomId
{
public:
  int	_Molecule;
  int	_Residue;
  int	_Atom;
public:
  static const int UndefId = -1;
public:
  AtomId() :
      _Molecule(UndefId),
      _Residue(UndefId),
      _Atom(UndefId) {};

  AtomId(uint moleculeId, uint residueId, chem::ConstitutionAtomIndex0N atomIndex) :
      _Molecule(moleculeId), _Residue(residueId), _Atom(atomIndex) {};

  inline bool operator==(AtomId const& other) const
  {
    if ( this->_Molecule != other._Molecule ) return false;
    if ( this->_Residue != other._Residue ) return false;
    return ( this->_Atom == other._Atom );
  }

  inline bool operator!=(AtomId const& other) const
  {
    if ( this->_Molecule != other._Molecule ) return true;
    if ( this->_Residue != other._Residue ) return true;
    return ( this->_Atom != other._Atom );
  }

  inline bool operator<(AtomId const& other) const
  {
    if ( this->_Molecule < other._Molecule ) return true;
    if ( this->_Residue < other._Residue ) return true;
    return ( this->_Atom < other._Atom );
  }




  bool isDefined() const { return this->_Molecule != UndefId; };
  int moleculeId() const { return this->_Molecule;};
  int residueId() const { return this->_Residue;};
  int atomId() const { return this->_Atom;};

  virtual ~AtomId() {};

  string asString() const;
};

};


namespace translate
{
template <>
struct	from_object<const chem::AtomId&>
{
  typedef	const chem::AtomId&		ExpectedType;
  typedef	chem::AtomId			DeclareType;
  DeclareType _v;
  from_object(core::T_sp o) : _v(core::clasp_to_fixnum(core::oCar(o)),
                                 core::clasp_to_fixnum(core::oCadr(o)),
                                 core::clasp_to_fixnum(core::oCaddr(o)))
  {
  }
};

template <>
struct	from_object<chem::AtomId>
{
  typedef	chem::AtomId		ExpectedType;
  typedef	chem::AtomId			DeclareType;
  DeclareType _v;
  from_object(core::T_sp o) : _v(core::clasp_to_fixnum(core::oCar(o)),
                                 core::clasp_to_fixnum(core::oCadr(o)),
                                 core::clasp_to_fixnum(core::oCaddr(o)))
  {
  }
};




template <>
struct	to_object<chem::AtomId>
{
  typedef	core::List_sp		ExpectedType;
  typedef	core::List_sp		DeclareType;
  static core::T_sp convert(chem::AtomId v)
  {_G();
    ExpectedType res = core::Cons_O::createList(
        core::clasp_make_fixnum(v.moleculeId()),
        core::clasp_make_fixnum(v.residueId()),
        core::clasp_make_fixnum(v.atomId()));
    return (res);
  }
};
};

#endif // chem_atomId_H
