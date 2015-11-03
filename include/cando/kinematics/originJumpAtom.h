#ifndef	kinematics_originJumpAtom_H
#define kinematics_originJumpAtom_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/atom.h>
#include <cando/kinematics/jump.h>
#include <cando/kinematics/jumpAtom.h>
#include <cando/chem/atomId.h>

namespace kinematics
{

    class OriginJumpAtom : public JumpAtom
    {
    public:
	static const NodeType nodeType = originJumpAtom;
    protected:
    public:
	/*! Empty ctor */
	OriginJumpAtom() {};

    OriginJumpAtom(const chem::AtomId& atomId,const string& comment) : JumpAtom(atomId,comment) {};



	virtual core::Symbol_sp typeSymbol() const;


	/*! The stubAtoms of an OriginJumpAtom are all myself */
	RefCountedAtomHandle stubAtom1() const { return this->_Me;};
	RefCountedAtomHandle stubAtom2() const { return this->_Me;};
	RefCountedAtomHandle stubAtom3(AtomTree_sp at) const { return this->_Me;};

	/*! Update the internal coordinates */
	virtual void updateInternalCoords(Stub& stub,
					  bool const recursive,
					  AtomTree_sp at);




    };



};
#endif
