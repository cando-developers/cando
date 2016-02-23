/*
    File: openmmState.h
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
