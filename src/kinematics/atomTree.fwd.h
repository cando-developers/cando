#ifndef kinematics_atomTree_fwd_H
#define kinematics_atomTree_fwd_H




namespace kinematics
{



    typedef enum { unused, jumpAtom, originJumpAtom, bondedAtom, rootBondedAtom, delayedBondedAtom } NodeType;





    class AtomHolder;
    FORWARD(AtomTree);
    FORWARD(Atom);


    class AtomTreeWalkFunctor
    {
    protected:
	core::Lisp_sp	_Lisp;
    public:
    AtomTreeWalkFunctor(core::Lisp_sp lisp) : _Lisp(lisp) {};

	core::Lisp_sp lisp() const { return this->_Lisp;};
	virtual void operator()(kinematics::Atom_sp atom) const;
    };


};
#endif
