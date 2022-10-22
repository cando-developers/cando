/*
    File: openmmState.cc
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
#define	DEBUG_LEVEL_NONE

#include <foundation.h>
#include <object.h>
#include <lisp.h>
#include <cando/openmmPackage/openmmFoundation.h>
#include <coordinateArray.h>
#include <clasp/core/numbers.h>
#include <cando/chem/atom.h>
#include <cando/chem/loop.h>
#include <cando/chem/forceField.h>
#include <cando/chem/bond.h>
#include <cando/chem/angle.h>
#include <cando/chem/properTorsion.h>
#include <cando/chem/improperTorsion.h>
#include <cando/units/unitsPackage.h>
#include <cando/units/quantity.h>
#include <cando/openmmPackage/openmmExpose.h>
#include <cando/openmmPackage/openmmState.h>
#include <wrappers.h>



namespace omm
{






    units::Quantity_sp State_O::getTime() const
    {
	IMPLEMENT_ME();
//	return this->_State.getTime();
    };

    void State_O::updateCandoMatterWithPositions(chem::Matter_sp matter) const
    {
	units::Quantity_sp coords;
	try
	{
	    coords = this->getPositions();
	}
	catch (...)
	{
	    SIMPLE_ERROR(("Could not get positions from openmm::state"));
	}
	chem::Loop lAtoms(matter,ATOMS);
	core::Symbol_sp particleIndexSymbol = _lisp->internWithPackageName(Pkg(),ParticleIndex);
	while ( lAtoms.advance() )
	{
	    chem::Atom_sp atom = lAtoms.getAtom();
	    int index = atom->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
	    Vector3 pos = coords->getElement_in_unit(index,_lisp->symbol(units::_sym_angstroms)->dyn<units::Unit_O>()).as<core::OVector3_O>()->get();
	    atom->setPosition(pos);
	}
    }

    units::Quantity_sp State_O::getPositions() const
    {
	core::SimpleVectorCoordinate_sp pos_in_nm; 
	try
	{
	    pos_in_nm = translate::to_object<const std::vector< OpenMM::Vec3 >& >::convert(_lisp,this->_State.getPositions());;
	}
	catch (...)
	{
	    SIMPLE_ERROR(("Could not get positions from openmm::state"));
	}
	units::Quantity_sp res =
	    units::Quantity_O::create(pos_in_nm,
				      _lisp->symbol(units::_sym_nanometers)->dyn<units::Unit_O>(),
				      _lisp);
        return res;
    };

    units::Quantity_sp State_O::getVelocities() const
    {
	core::SimpleVectorCoordinate_sp vel_in_nm_per_ps;
	try
	{
	    vel_in_nm_per_ps = translate::to_object<const std::vector< OpenMM::Vec3 >& >::convert(_lisp,this->_State.getVelocities());
	}
	catch (...)
	{
	    SIMPLE_ERROR(("Could not get velocities from openmm::state"));
	}
	units::Quantity_sp res = units::Quantity_O::create(vel_in_nm_per_ps,
							   _lisp->symbol(units::_sym_nm_per_ps)->dyn<units::Unit_O>(),
							   _lisp);
        return res;
    };

#if 0
    core::SimpleVectorCoordinate_sp State_O::getForces() const
    {
	// Get the forces in whatever units they come and convert them to SI 
	core::SimpleVectorCoordinate_sp res;
	try
	{
	    res = translate::to_object<const std::vector< OpenMM::Vec3 >& >::convert(_lisp,this->_State.getForces());
	}
	catch (...)
	{
	    SIMPLE_ERROR(("Could not get forces from openmm::state"));
	}
        return res;
    };
#endif




    units::Quantity_sp State_O::getKineticEnergy()
    {
	double e_in_kj_per_mole;
	try
	{
	    e_in_kj_per_mole = this->_State.getKineticEnergy();
	}
	catch (...)
	{
	    SIMPLE_ERROR(("Could not get kinetic energy from openmm::state"));
	}
	units::Quantity_sp f =
	    units::Quantity_O::create(e_in_kj_per_mole,
				      _lisp->symbol(units::_sym_kilojoules_per_mole)->dyn<units::Unit_O>(),
				      _lisp);
	return f;
    }

    units::Quantity_sp State_O::getPotentialEnergy()
    {
	double e_in_kj_per_mole;
	try
	{
	     e_in_kj_per_mole = this->_State.getPotentialEnergy();
	}
	catch (...)
	{
	    SIMPLE_ERROR(("Could not get potential energy from openmm::state"));
	}
	units::Quantity_sp f =
	    units::Quantity_O::create(e_in_kj_per_mole,
				      _lisp->symbol(units::_sym_kilojoules_per_mole)->dyn<units::Unit_O>(),
				      _lisp);
	ASSERTNOTNULL(f);
	return f;
    }

#if 0
    core::Cons_sp State_O::getPeriodicBoxVectors()
{
    OpenMM::Vec3 a,b,c;
    this->_State.getPeriodicBoxVectors(a,b,c);
    core::Cons_sp res = core::Cons_O::createList(translate::to_object<OpenMM::Vec3>::convert(_lisp,a),
						 translate::to_object<OpenMM::Vec3>::convert(_lisp,b),
						 translate::to_object<OpenMM::Vec3>::convert(_lisp,c),
				   _lisp);
    return res;
}
#endif




};
