#ifndef	kinematics_monomerId_H
#define kinematics_monomerId_H

#include <clasp/core/common.h>
#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>


namespace kinematics
{


    class MonomerId
    {
	friend class MonomerId_O;
    protected:
	int	_Chain;
	int	_Monomer;
    public:
	static const int UndefId = -1;
    public:
    MonomerId() :
	_Chain(UndefId),
	    _Monomer(UndefId) {};

    MonomerId(uint chainId, uint monomerId) :
	_Chain(chainId), _Monomer(monomerId) {};

	int chainId() const { return this->_Chain;};
	int monomerId() const { return this->_Monomer;};

	virtual ~MonomerId() {};

	string asString() const;
    };


    FORWARD(MonomerId);
    /*! @class Wrap MonomerId inside of a Lisp class */
    class MonomerId_O : public core::T_O
    {
	friend class FoldTree_O;
	LISP_BASE1(core::T_O);
	LISP_CLASS(kinematics,KinPkg,MonomerId_O,"MonomerId");
#if INIT_TO_FACTORIES
    public:
	static MonomerId_sp make(const int chain, const int monomer);
#else
	DECLARE_INIT();
#endif
//	DECLARE_ARCHIVE();
    protected:
	MonomerId	_MonomerId;
    public:
    MonomerId_O() : _MonomerId() {}
	virtual ~MonomerId_O() {};
	string asString() const {return this->_MonomerId.asString();};
    }; // MonomerId_O


};


namespace translate
{
    template <>
    struct	from_object<kinematics::MonomerId const&>
    {
	typedef	kinematics::MonomerId const&		ExpectedType;
	typedef	kinematics::MonomerId			DeclareType;
	DeclareType _v;
    from_object(core::T_sp o) : _v(core::oCar(o).as<core::Fixnum_I>().unsafe_fixnum(),
                                   core::oCadr(o).as<core::Fixnum_I>().unsafe_fixnum()) {}
    };


};



#endif // kinematics_monomerId_H
