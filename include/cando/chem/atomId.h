#ifndef	chem_atomId_H
#define chem_atomId_H

#include <clasp/core/common.h>
#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/chem/chemPackage.h>


namespace chem
{


    class AtomId
    {
	friend class AtomId_O;
    protected:
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


    FORWARD(AtomId);
    /*! @class Wrap AtomId inside of a Lisp class */
    class AtomId_O : public core::T_O
    {
	friend class AtomIdToAtomMap_O;
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,AtomId_O,"AtomId");
//	DECLARE_ARCHIVE();
    protected:
	AtomId	_AtomId;
    public:
	static AtomId_sp make(int mol, int res, int atom);
    public:
	int moleculeId() const { return this->_AtomId.moleculeId();};
	int residueId() const { return this->_AtomId.residueId();};
	int atomId() const { return this->_AtomId.atomId();};

	virtual ~AtomId_O() {};
	string asString() const {return this->_AtomId.asString();};
    }; // AtomId_O



};


namespace translate
{
    template <>
    struct	from_object<const chem::AtomId&>
    {
	typedef	const chem::AtomId&		ExpectedType;
	typedef	chem::AtomId			DeclareType;
	DeclareType _v;
    from_object(core::T_sp o) : _v(core::clasp_to_fixnum(gc::As<core::Fixnum_sp>(core::oCar(o))),
                                   core::clasp_to_fixnum(core::oCadr(o).as<core::Fixnum_I>()),
                                   core::clasp_to_fixnum(core::oCaddr(o).as<core::Fixnum_I>()))
	{
#if 0
	    core::List_sp args = o.as<core::Cons_O>();
	    int moleculeId = core::oCar(args).as<core::Rational_O>()->as_int();
	    int residueId = args->ocadr().as<core::Rational_O>()->as_int();
	    int atomId = args->ocaddr().as<core::Rational_O>()->as_int();
	    _v.set(moleculeId,residueId,atomId);
#endif
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

TRANSLATE(chem::AtomId_O);


#endif // chem_atomId_H
