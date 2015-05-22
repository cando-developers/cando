#ifndef	openmmState_H //[
#define openmmState_H

#include <OpenMM.h>


#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <foundation.h>
#include <object.h>
#include <ovector3.h>
#include <cando/units/quantity.fwd.h>
#include <cando/chem/matter.h>
#include <cando/openmmPackage/openmmFoundation.h>
#include <externalObject.h>



namespace omm
{


    class Context_O;

    class State_O : public core::T_O
    {
	friend class Context_O;
	LISP_BASE1(core::T_O);
	LISP_CLASS(OmmPkg,State_O,"State");
    protected:
	OpenMM::State	_State;
    public:
	units::Quantity_sp getTime() const;
	/*! Fill the atom positions (in angstroms) from the State */
	void updateCandoMatterWithPositions(chem::Matter_sp matter) const;

	/*! Return a positions array with units */
	units::Quantity_sp getPositions() const;

	/*! Return a velocities array with units */
	units::Quantity_sp getVelocities() const;
	units::Quantity_sp getForces() const;
	units::Quantity_sp getKineticEnergy();
	units::Quantity_sp getPotentialEnergy();
	core::Cons_sp getPeriodicBoxVectors();


	DEFAULT_CTOR_DTOR(State_O);
    };

};

TRANSLATE(omm::State_O);

#endif //]
