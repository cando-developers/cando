/*
    File: monomerId.h
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
    class MonomerId_O : public core::General_O
    {
	friend class FoldTree_O;
	LISP_CLASS(kinematics,KinPkg,MonomerId_O,"MonomerId",core::General_O);
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
