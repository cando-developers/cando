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
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,ZMatrixInternal_O,"ZMatrixInternal");
    public:
	void initialize();
//	void	archiveBase(core::ArchiveP node);
    private:
	ZMatrix_wp	_WeakZMatrix;
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
	virtual double	getValue() { return this->_Value; };

	ZMatrixInternal_O( const ZMatrixInternal_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ZMatrixInternal_O);
    };





    SMART(ZMatrixBondInternal);
    class ZMatrixBondInternal_O : public ZMatrixInternal_O
    {
	LISP_BASE1(ZMatrixInternal_O);
	LISP_CLASS(chem,ChemPkg,ZMatrixBondInternal_O,"ZMatrixBondInternal");
    public:
//	void	archiveBase(core::ArchiveP node);
    private:
	uint	_AtomBond;
    public:
	static ZMatrixBondInternal_sp create(Atom_sp newAtom, Atom_sp bondToAtom,
					     ZMatrixInternal_O::atomMap atomIndices, ZMatrix_sp zMatrix );
    public:
	void setup(Atom_sp newAtom, Atom_sp bondToAtom,
		   ZMatrixInternal_O::atomMap atomIndices);
	string	getBondAtomZMatrixName();
	int	getBondAtomIndex() const { return this->_AtomBond;};

	virtual string	__repr__() const;

	virtual	void	extractInternal();
	virtual	void	constrainIfNoAtomsWithMask(uint m);

	ZMatrixBondInternal_O( const ZMatrixBondInternal_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ZMatrixBondInternal_O);
    };

    inline ZMatrixBondInternal_sp ZMatrixBondInternal_O::create(Atom_sp newAtom, Atom_sp bondToAtom,
								ZMatrixInternal_O::atomMap atomIndices, ZMatrix_sp zMatrix )
    {_G();
	GC_ALLOCATE(ZMatrixBondInternal_O, zz ); // RP_Create<ZMatrixBondInternal_O>(e);
	zz->setZMatrix(zMatrix);
	zz->setup(newAtom,bondToAtom,atomIndices);
	return zz;
    };



    SMART(ZMatrixAngleInternal);
    class ZMatrixAngleInternal_O : public ZMatrixInternal_O
    {
	LISP_BASE1(ZMatrixInternal_O);
	LISP_CLASS(chem,ChemPkg,ZMatrixAngleInternal_O,"ZMatrixAngleInternal");

    public:
//	void	archiveBase(core::ArchiveP node);
    private:
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
	string	getBondAtomZMatrixName();
	string	getAngleAtomZMatrixName();
	int getBondAtomIndex() const { return this->_AtomBond;}
	int getAngleAtomIndex() const { return this->_AtomAngle;}

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
	GC_ALLOCATE(ZMatrixAngleInternal_O, zz ); // RP_Create<ZMatrixAngleInternal_O>(e);
	zz->setZMatrix(zMatrix);
	zz->setup(newAtom,bondToAtom,angleToAtom,atomIndices);
	return zz;
    };


    SMART(ZMatrixDihedralInternal);
    class ZMatrixDihedralInternal_O : public ZMatrixInternal_O
    {
	LISP_BASE1(ZMatrixInternal_O);
	LISP_CLASS(chem,ChemPkg,ZMatrixDihedralInternal_O,"ZMatrixDihedralInternal");

    public:
//	void	archiveBase(core::ArchiveP node);
    private:
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
	string	getBondAtomZMatrixName();
	string	getAngleAtomZMatrixName();
	string	getDihedralAtomZMatrixName();
	int getBondAtomIndex() const { return this->_AtomBond;}
	int getAngleAtomIndex() const { return this->_AtomAngle;}
	int getDihedralAtomIndex() const { return this->_AtomDihedral;}

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
	GC_ALLOCATE(ZMatrixDihedralInternal_O, zz ); // RP_Create<ZMatrixDihedralInternal_O>(e);
	zz->setZMatrix(zMatrix);
	zz->setup(newAtom,bondToAtom,angleToAtom,dihedralToAtom,atomIndices);
	return zz;
    };




    SMART(ZMatrixEntry);
    class ZMatrixEntry_O : public core::CxxObject_O
    {
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,ZMatrixEntry_O,"ZMatrixEntry");
    public:
	void initialize();
	friend class ZMatrix_O;
    public:
//	void	archiveBase(core::ArchiveP node);
    protected:
	string				_ZMatrixAtomName;
	Atom_sp				_Atom;
	ZMatrixBondInternal_sp		_Bond;
	ZMatrixAngleInternal_sp		_Angle;
	ZMatrixDihedralInternal_sp	_Dihedral;
    public:
	static ZMatrixEntry_sp create(Atom_sp atom, ZMatrixInternal_O::atomMap atomIndices);
	void setup(Atom_sp atom, ZMatrixInternal_O::atomMap atomIndices);
    public:
	Atom_sp	getAtom() { return this->_Atom;};
	string	getTargetAtomZMatrixName() { return this->_ZMatrixAtomName; };

	bool	hasBond() {_OF(); ASSERTNOTNULL(this->_Bond);return this->_Bond.notnilp(); };
	ZMatrixBondInternal_sp	getBond() {_OF(); ASSERTNOTNULL(this->_Bond);return this->_Bond;};
	bool	hasAngle() {_OF(); ASSERTNOTNULL(this->_Angle);return this->_Angle.notnilp();};
	ZMatrixAngleInternal_sp	getAngle() {_OF(); ASSERTNOTNULL(this->_Angle);return this->_Angle;};
	bool	hasDihedral() {_OF(); ASSERTNOTNULL(this->_Dihedral);return this->_Dihedral.notnilp(); };
	ZMatrixDihedralInternal_sp	getDihedral() {_OF(); ASSERTNOTNULL(this->_Dihedral);return this->_Dihedral;};

	ZMatrixEntry_O( const ZMatrixEntry_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ZMatrixEntry_O);
    };



    SMART(ZMatrix);
    class ZMatrix_O : public core::CxxObject_O
    {
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,ZMatrix_O,"ZMatrix");
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
//	void	archiveBase(core::ArchiveP node);
    private:
	Matter_sp		_Matter;
        gctools::Vec0<ZMatrixEntry_sp>	_ZMatrix;
        gctools::Vec0<ZMatrixInternal_sp>	_Internals;
    public:
	typedef	gctools::Vec0<ZMatrixEntry_sp>::iterator		zMatrixEntryIterator;
	typedef	gctools::Vec0<ZMatrixEntry_sp>::const_iterator		const_zMatrixEntryIterator;
	typedef	gctools::Vec0<ZMatrixInternal_sp>::iterator	internalIterator;
    protected:
	Atom_sp	_getAtomAtIndex(uint i);
	string	_getAtomZMatrixNameAtIndex(uint i) const;
	void	defineForMatterWithStartingAtom(Matter_sp matter, Atom_sp atom);

    public:
	void	defineForMatter(Matter_sp matter);

	/*! Extract the internal coordinates for the ZMatrix */
	void	extractInternals();

	/*! Get the atom name for a particular entry */
	string getZMatrixAtomNameAtIndex(uint i) const;

	/*! Return the entries as a Cons */
	core::List_sp entriesAsCons() const;

	const_zMatrixEntryIterator	begin_ZMatrixEntries() const	{ return this->_ZMatrix.begin(); };
	const_zMatrixEntryIterator	end_ZMatrixEntries() const	{ return this->_ZMatrix.end(); };
	zMatrixEntryIterator	begin_ZMatrixEntries() 	{ return this->_ZMatrix.begin(); };
	zMatrixEntryIterator	end_ZMatrixEntries() 	{ return this->_ZMatrix.end(); };

	internalIterator	begin_Internals()	{ return this->_Internals.begin(); };
	internalIterator	end_Internals()	{ return this->_Internals.end(); };

	void	constrainInternalsWithNoAtomsWithMask(uint mask);


	ZMatrixEntry_sp zMatrixEntryWithName(string const& atomZMatrixName) const;

	string allZMatrixAtomNamesAsString() const;

	ZMatrix_O( const ZMatrix_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(ZMatrix_O);
    };




};

    TRANSLATE(chem::ZMatrixInternal_O);
    TRANSLATE(chem::ZMatrixBondInternal_O);
    TRANSLATE(chem::ZMatrixAngleInternal_O);
    TRANSLATE(chem::ZMatrixDihedralInternal_O);
    TRANSLATE(chem::ZMatrixEntry_O);
    TRANSLATE(chem::ZMatrix_O);



#endif //]
