/*
    File: cipPrioritizer.h
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
#ifndef	CipPrioritizer_H //[
#define CipPrioritizer_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>


namespace chem {

SMART(Matter);
SMART(Atom);



struct AtomPriority {
  AtomPriority(Atom_sp a, int p) : _Atom(a), _relativePriority(p) {};
  AtomPriority() : _Atom(nil<Atom_O>()), _relativePriority(0) {};
  Atom_sp _Atom;
  int     _relativePriority;
};


SMART(CipPrioritizer );
class CipPrioritizer_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,CipPrioritizer_O,"CipPrioritizer",core::CxxObject_O);
private:
	// instance variables
  gctools::Vec0<int>		        _p;
  gctools::Vec0< gctools::Vec0<int> > 	_s;
public:
  CipPrioritizer_O() {};
  
	void	initialize();
public:
    static core::HashTable_sp assignPrioritiesHashTable(Matter_sp molOrAgg);
  CL_DEF_CLASS_METHOD static CipPrioritizer_sp make_CipPrioritizer() {
    CipPrioritizer_sp cip = gc::GC<CipPrioritizer_O>::allocate();
    return cip;
  }
public:
  int	getP(AtomPriority& a);
  gctools::Vec0<int>& getS(AtomPriority& a);
  void assignCahnIngoldPrelogPriorityToAtomsRelativePriority(Matter_sp molOrAgg,core::HashTable_sp cip);
public:

  core::HashTable_sp setStereochemicalTypeForAllAtoms(Matter_sp molOrAgg);
  CipPrioritizer_O( const CipPrioritizer_O& ss ); //!< Copy constructor

};


};
#endif //]
