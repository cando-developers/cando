/*
    File: zMatrix.h
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
#ifndef	ZMatrix_H //[
#define ZMatrix_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/adapt/stringList.h>
#include <cando/chem/chemPackage.h>
//#include "conformationCollection.h"





namespace chem
{

  SMART(Atom);
    SMART(ZMatrix);

    SMART(ZMatrixInternal );
    class ZMatrixInternal_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,ZMatrixInternal_O,"ZMatrixInternal",core::CxxObject_O);
    public:
	void initialize();
    private:
	ZMatrix_sp	_WeakZMatrix;
	string		_InternalName;
	bool		_Constrain;
	uint	_AtomNew;
	double		_Value;
    public:
	typedef gctools::SmallMap<Atom_sp,uint> atomMap; // m a p<Atom_sp,uint>	atomMap;

    public:
	void	setZMatrix(ZMatrix_sp zm) { this->_WeakZMatrix = zm; };
	ZMatrix_sp	getZMatrix() {_OF(); ASSERTNOTNULL(this->_WeakZMatrix);return this->_WeakZMatrix;};

	void	setAtomNew(Atom_sp atom, atomMap atomIndices);
	void	setInternalName(const string& name) { this->_InternalName = name; };
	bool	getConstrain() { return this->_Constrain;};
	void	setConstrain(bool c) { this->_Constrain = c;};
	virtual	void	extractInternal() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
	virtual	void	constrainIfNoAtomsWithMask(uint m) {_OF(); SUBCLASS_MUST_IMPLEMENT();};

	virtual string __repr__() const;

	string	getInternalName() { return this->_InternalName;};
	Atom_sp	getAtomNew();

	void	setValue(double v) { this->_Value = v;};
CL_LISPIFY_NAME("getValue");
CL_DEFMETHOD 	virtual double	getValue() { return this->_Value; };

	ZMatrixInternal_O( const ZMatrixInternal_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ZMatrixInternal_O);
    };





    SMART(ZMatrixBondInternal);
    class ZMatrixBondInternal_O : public ZMatrixInternal_O
    {
	LISP_CLASS(chem,ChemPkg,ZMatrixBondInternal_O,"ZMatrixBondInternal",ZMatrixInternal_O);
    public:
	uint	_AtomBond;
    public:
	static ZMatrixBondInternal_sp create(Atom_sp newAtom, Atom_sp bondToAtom,
					     ZMatrixInternal_O::atomMap atomIndices, ZMatrix_sp zMatrix );
    public:
	void setup(Atom_sp newAtom, Atom_sp bondToAtom,
		   ZMatrixInternal_O::atomMap atomIndices);
        core::T_sp	getBondAtomZMatrixName();
CL_LISPIFY_NAME("getBondAtomIndex");
CL_DEFMETHOD 	int	getBondAtomIndex() const { return this->_AtomBond;};

	virtual string	__repr__() const;

	virtual	void	extractInternal();
	virtual	void	constrainIfNoAtomsWithMask(uint m);

	ZMatrixBondInternal_O( const ZMatrixBondInternal_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ZMatrixBondInternal_O);
    };

    inline ZMatrixBondInternal_sp ZMatrixBondInternal_O::create(Atom_sp newAtom, Atom_sp bondToAtom,
								ZMatrixInternal_O::atomMap atomIndices, ZMatrix_sp zMatrix )
    {_G();
      auto  zz  = gctools::GC<ZMatrixBondInternal_O>::allocate_with_default_constructor();
	zz->setZMatrix(zMatrix);
	zz->setup(newAtom,bondToAtom,atomIndices);
	return zz;
    };



    SMART(ZMatrixAngleInternal);
    class ZMatrixAngleInternal_O : public ZMatrixInternal_O
    {
	LISP_CLASS(chem,ChemPkg,ZMatrixAngleInternal_O,"ZMatrixAngleInternal",ZMatrixInternal_O);

    public:
	uint	_AtomBond;
	uint	_AtomAngle;
    public:
	static ZMatrixAngleInternal_sp create(Atom_sp newAtom, Atom_sp bondToAtom,
					      Atom_sp angleToAtom,
					      ZMatrixInternal_O::atomMap atomIndices,
					      ZMatrix_sp zMatrix );
    public:
    	void setup(Atom_sp newAtom, Atom_sp bondToAtom,
		   Atom_sp angleToAtom,
		   ZMatrixInternal_O::atomMap atomIndices);
    public:
        core::T_sp	getBondAtomZMatrixName();
        core::T_sp	getAngleAtomZMatrixName();
CL_LISPIFY_NAME("getBondAtomIndex");
CL_DEFMETHOD 	int getBondAtomIndex() const { return this->_AtomBond;}
CL_LISPIFY_NAME("getAngleAtomIndex");
CL_DEFMETHOD 	int getAngleAtomIndex() const { return this->_AtomAngle;}

	virtual string	__repr__() const;

	virtual	void	extractInternal();
	virtual	void	constrainIfNoAtomsWithMask(uint m);

	virtual double	getValue() { return this->ZMatrixInternal_O::getValue()/0.0174533;};

	ZMatrixAngleInternal_O( const ZMatrixAngleInternal_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ZMatrixAngleInternal_O);
    };

    inline ZMatrixAngleInternal_sp ZMatrixAngleInternal_O::create(Atom_sp newAtom, Atom_sp bondToAtom,
								  Atom_sp angleToAtom,
								  ZMatrixInternal_O::atomMap atomIndices, ZMatrix_sp zMatrix)
    {_G();
      auto  zz  = gctools::GC<ZMatrixAngleInternal_O>::allocate_with_default_constructor();
	zz->setZMatrix(zMatrix);
	zz->setup(newAtom,bondToAtom,angleToAtom,atomIndices);
	return zz;
    };


    SMART(ZMatrixDihedralInternal);
    class ZMatrixDihedralInternal_O : public ZMatrixInternal_O
    {
	LISP_CLASS(chem,ChemPkg,ZMatrixDihedralInternal_O,"ZMatrixDihedralInternal",ZMatrixInternal_O);

    public:
	uint	_AtomBond;
	uint	_AtomAngle;
	uint	_AtomDihedral;
    public:
	static ZMatrixDihedralInternal_sp create(Atom_sp newAtom, Atom_sp bondToAtom,
						 Atom_sp angleToAtom, Atom_sp dihedralToAtom,
						 ZMatrixInternal_O::atomMap atomIndices, ZMatrix_sp zMatrix );
						 
    public:
    	void setup(Atom_sp newAtom, Atom_sp bondToAtom,
		   Atom_sp angleToAtom, Atom_sp dihedralToAtom,
		   ZMatrixInternal_O::atomMap atomIndices);

    public:
        core::T_sp	getBondAtomZMatrixName();
        core::T_sp	getAngleAtomZMatrixName();
        core::T_sp	getDihedralAtomZMatrixName();
CL_LISPIFY_NAME("getBondAtomIndex");
CL_DEFMETHOD 	int getBondAtomIndex() const { return this->_AtomBond;}
CL_LISPIFY_NAME("getAngleAtomIndex");
CL_DEFMETHOD 	int getAngleAtomIndex() const { return this->_AtomAngle;}
CL_LISPIFY_NAME("getDihedralAtomIndex");
CL_DEFMETHOD 	int getDihedralAtomIndex() const { return this->_AtomDihedral;}

	virtual string	__repr__() const;

	virtual	void	extractInternal();
	virtual	void	constrainIfNoAtomsWithMask(uint m);
	virtual double	getValue() { return this->ZMatrixInternal_O::getValue()/0.0174533;};

	ZMatrixDihedralInternal_O( const ZMatrixDihedralInternal_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ZMatrixDihedralInternal_O);
    };

    inline ZMatrixDihedralInternal_sp ZMatrixDihedralInternal_O::create(Atom_sp newAtom, Atom_sp bondToAtom,
									Atom_sp angleToAtom, Atom_sp dihedralToAtom,
									ZMatrixInternal_O::atomMap atomIndices, ZMatrix_sp zMatrix)
    {_G();
      auto  zz  = gctools::GC<ZMatrixDihedralInternal_O>::allocate_with_default_constructor();
	zz->setZMatrix(zMatrix);
	zz->setup(newAtom,bondToAtom,angleToAtom,dihedralToAtom,atomIndices);
	return zz;
    };




    SMART(ZMatrixEntry);
    class ZMatrixEntry_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,ZMatrixEntry_O,"ZMatrixEntry",core::CxxObject_O);
    public:
	void initialize();
	friend class ZMatrix_O;
    public:
      core::T_sp         _ZMatrixAtomName;
      Atom_sp				_Atom;
      ZMatrixBondInternal_sp		_Bond;
      ZMatrixAngleInternal_sp		_Angle;
      ZMatrixDihedralInternal_sp	_Dihedral;
    public:
	static ZMatrixEntry_sp create(Atom_sp atom, ZMatrixInternal_O::atomMap atomIndices);
	void setup(Atom_sp atom, ZMatrixInternal_O::atomMap atomIndices);
    public:
CL_LISPIFY_NAME("getAtom");
CL_DEFMETHOD 	Atom_sp	getAtom() { return this->_Atom;};
CL_LISPIFY_NAME("getTargetAtomZMatrixName");
 CL_DEFMETHOD 	core::T_sp	getTargetAtomZMatrixName() { return this->_ZMatrixAtomName; };

CL_LISPIFY_NAME("hasBond");
CL_DEFMETHOD 	bool	hasBond() {_OF(); ASSERTNOTNULL(this->_Bond);return this->_Bond.notnilp(); };
CL_LISPIFY_NAME("getBond");
CL_DEFMETHOD 	ZMatrixBondInternal_sp	getBond() {_OF(); ASSERTNOTNULL(this->_Bond);return this->_Bond;};
CL_LISPIFY_NAME("hasAngle");
CL_DEFMETHOD 	bool	hasAngle() {_OF(); ASSERTNOTNULL(this->_Angle);return this->_Angle.notnilp();};
CL_LISPIFY_NAME("getAngle");
CL_DEFMETHOD 	ZMatrixAngleInternal_sp	getAngle() {_OF(); ASSERTNOTNULL(this->_Angle);return this->_Angle;};
CL_LISPIFY_NAME("hasDihedral");
CL_DEFMETHOD 	bool	hasDihedral() {_OF(); ASSERTNOTNULL(this->_Dihedral);return this->_Dihedral.notnilp(); };
CL_LISPIFY_NAME("getDihedral");
CL_DEFMETHOD 	ZMatrixDihedralInternal_sp	getDihedral() {_OF(); ASSERTNOTNULL(this->_Dihedral);return this->_Dihedral;};

	ZMatrixEntry_O( const ZMatrixEntry_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ZMatrixEntry_O);
    };



    SMART(ZMatrix);
    class ZMatrix_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,ZMatrix_O,"ZMatrix",core::CxxObject_O);
#if INIT_TO_FACTORIES
    public:
	static ZMatrix_sp make(Matter_sp matter, Atom_sp roomAtom);
#else
	DECLARE_INIT();
#endif
	friend class	ZMatrixInternal_O;
	friend class	ZMatrixBondInternal_O;
	friend class	ZMatrixAngleInternal_O;
	friend class	ZMatrixDihedralInternal_O;
    public:
	Matter_sp		_Matter;
        gctools::Vec0<ZMatrixEntry_sp>	_ZMatrix;
        gctools::Vec0<ZMatrixInternal_sp>	_Internals;
    public:
	typedef	gctools::Vec0<ZMatrixEntry_sp>::iterator		zMatrixEntryIterator;
	typedef	gctools::Vec0<ZMatrixEntry_sp>::const_iterator		const_zMatrixEntryIterator;
	typedef	gctools::Vec0<ZMatrixInternal_sp>::iterator	internalIterator;
    protected:
	Atom_sp	_getAtomAtIndex(uint i);
        core::T_sp	_getAtomZMatrixNameAtIndex(uint i) const;
	void	defineForMatterWithStartingAtom(Matter_sp matter, Atom_sp atom);

    public:
	void	defineForMatter(Matter_sp matter);

	/*! Extract the internal coordinates for the ZMatrix */
	void	extractInternals();

	/*! Get the atom name for a particular entry */
        core::T_sp getZMatrixAtomNameAtIndex(uint i) const;

	/*! Return the entries as a Cons */
	core::List_sp entriesAsList() const;

	const_zMatrixEntryIterator	begin_ZMatrixEntries() const	{ return this->_ZMatrix.begin(); };
	const_zMatrixEntryIterator	end_ZMatrixEntries() const	{ return this->_ZMatrix.end(); };
	zMatrixEntryIterator	begin_ZMatrixEntries() 	{ return this->_ZMatrix.begin(); };
	zMatrixEntryIterator	end_ZMatrixEntries() 	{ return this->_ZMatrix.end(); };

	internalIterator	begin_Internals()	{ return this->_Internals.begin(); };
	internalIterator	end_Internals()	{ return this->_Internals.end(); };

	void	constrainInternalsWithNoAtomsWithMask(uint mask);


	ZMatrixEntry_sp zMatrixEntryWithName(core::T_sp atomZMatrixName) const;

	string allZMatrixAtomNamesAsString() const;

	ZMatrix_O( const ZMatrix_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(ZMatrix_O);
    };




};




#endif //]
