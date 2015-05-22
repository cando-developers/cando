#ifndef	kinematics_rootBondedAtom_H
#define kinematics_rootBondedAtom_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/kinFoundation.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/rootAtomInfo.h>
#include <cando/kinematics/bondedAtom.h>



namespace kinematics
{

    class RootBondedAtom : public BondedAtom
    {
    public:
	static const NodeType nodeType = rootBondedAtom;
    protected:
	/*! Store the Id of the Bond1 atom if we have one
	  otherwise it will be UndefinedUnsignedInt */
	RootAtomInfo	_RootInfo;
    private:
    public:
    RootBondedAtom() : BondedAtom() {};
    RootBondedAtom(const chem::AtomId& atomId,const string& comment) :
	BondedAtom(atomId,comment) {};


	virtual RootAtomInfo const * rootAtomInfo() const { return &this->_RootInfo;};


	/*! Set the Bond1Id for this atom */

	void setup(core::Symbol_sp constitutionName,
		   core::Symbol_sp topologyName,
		   chem::Plug_sp inPlug);

	virtual core::Symbol_sp typeSymbol() const;

    };



};
#endif
