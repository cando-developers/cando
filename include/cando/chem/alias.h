/*
    File: alias.h
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
       
       
#ifndef	Alias_H //[
#define Alias_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>


#include <cando/chem/chemPackage.h>


namespace chem
{


    SMART(ScorerState);

    SMART(Alias);
    class Alias_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,Alias_O,"Alias",core::CxxObject_O);
#if INIT_TO_FACTORIES
    public:
	static Alias_sp make(core::Symbol_sp monomerSym, core::Symbol_sp atomSym);
#else
#endif
    public:
    public:
	void	initialize();
    public:
	static Alias_sp create(core::Symbol_sp monomerAlias, core::Symbol_sp atomAlias )
	{
	    GC_ALLOCATE(Alias_O, a ); // RP_Create<Alias_O>(e);
	    a->setMonomerAlias(monomerAlias);
	    a->setAtomAlias(atomAlias);
	    return a;
	}
	static core::T_sp prim_alias(core::Function_sp, core::List_sp args);
    private:
	core::Symbol_sp	_MonomerAlias;
	core::Symbol_sp	_AtomAlias;
    private: // Do not archive indices
	uint	_MonomerId;
	uint	_AtomId;
    public:

	core::Symbol_sp getMonomerAlias() { return this->_MonomerAlias;};
	void	setMonomerAlias(core::Symbol_sp a) { this->_MonomerAlias = a;};
	core::Symbol_sp getAtomAlias() { return this->_AtomAlias;};
	void	setAtomAlias(core::Symbol_sp a) { this->_AtomAlias = a;};


#if 0
	/*! If the alchemist oligomer changes then update this
	 * alias to the correct monomer index
	 */
	void oligomerChanged(PointProvider_sp alchemist);
	/*! If the alchemist oligomer changes then update this
	 * alias to the correct atom index
	 */
	void sequenceChanged(PointProvider_sp alchemist);
	/*! Get the atom position for this alias
	 */
//	Vector3 getAtomPosition(PointProvider_sp alchemist);

	/*! Get the transformed atom position for this alias
	 */
	Vector3& getUntransformedAtomPosition(PointProvider_sp alchemist);
#endif


	virtual bool equal(core::T_sp obj) const;

	string description() const;
	string __repr__() const { return this->_MonomerAlias->__repr__()+"@"+this->_AtomAlias->__repr__();};
	Alias_O(const core::Lisp_sp& lisp);
	Alias_O( const Alias_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(Alias_O);
    };




};
TRANSLATE(chem::Alias_O);
#endif //]
