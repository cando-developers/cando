#define	DEBUG_LEVEL_FULL

#include "foundation.h"
#include "object.h"
#include "lisp.h"
#include "openmmFoundation.h"
#include "coordinateArray.h"
#include "core/numbers.h"
#include "chem/atom.h"
#include "chem/loop.h"
#include "chem/forceField.h"
#include "chem/bond.h"
#include "chem/angle.h"
#include "chem/properTorsion.h"
#include "chem/improperTorsion.h"
#include "units/unitsPackage.h"
#include "units/quantity.h"
#include "openmmExpose.h"
#include "openmmState.h"
#include "wrappers.h"



namespace omm
{


    void State_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<State_O>()
	    .def("getTime",&State_O::getTime)
	    .def("getPositions",&State_O::getPositions)
	    .def("getVelocities",&State_O::getVelocities)
//	    .def("getForces",&State_O::getForces)
	    .def("getKineticEnergy",&State_O::getKineticEnergy)
	    .def("getPotentialEnergy",&State_O::getPotentialEnergy)
//	    .def("getPeriodicBoxVectors",&State_O::getPeriodicBoxVectors)
	    .def("updateCandoMatterWithPositions",&State_O::updateCandoMatterWithPositions)
	    ;
    }


    void State_O::exposePython(core::Lisp_sp lisp)
    {_G();
	PYTHON_CLASS(OmmPkg,State,"","",_lisp)
	    .def("getTime",&State_O::getTime)
	    .def("getPositions",&State_O::getPositions)
	    .def("getVelocities",&State_O::getVelocities)
//	    .def("getForces",&State_O::getForces)
	    .def("getKineticEnergy",&State_O::getKineticEnergy)
	    .def("getPotentialEnergy",&State_O::getPotentialEnergy)
//	    .def("getPeriodicBoxVectors",&State_O::getPeriodicBoxVectors)
	    .def("updateCandoMatterWithPositions",&State_O::updateCandoMatterWithPositions)
	    ;
    }


    units::Quantity_sp State_O::getTime() const
    {
	IMPLEMENT_ME();
//	return this->_State.getTime();
    };

    void State_O::updateCandoMatterWithPositions(chem::Matter_sp matter) const
    {_OF();
	units::Quantity_sp coords;
	try
	{
	    coords = this->getPositions();
	}
	catch (...)
	{
	    SIMPLE_ERROR(BF("Could not get positions from openmm::state"));
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
    {_OF();
	core::CoordinateArray_sp pos_in_nm; 
	try
	{
	    pos_in_nm = translate::to_object<const std::vector< OpenMM::Vec3 >& >::convert(_lisp,this->_State.getPositions());;
	}
	catch (...)
	{
	    SIMPLE_ERROR(BF("Could not get positions from openmm::state"));
	}
	units::Quantity_sp res =
	    units::Quantity_O::create(pos_in_nm,
				      _lisp->symbol(units::_sym_nanometers)->dyn<units::Unit_O>(),
				      _lisp);
        return res;
    };

    units::Quantity_sp State_O::getVelocities() const
    {_OF();
	core::CoordinateArray_sp vel_in_nm_per_ps;
	try
	{
	    vel_in_nm_per_ps = translate::to_object<const std::vector< OpenMM::Vec3 >& >::convert(_lisp,this->_State.getVelocities());
	}
	catch (...)
	{
	    SIMPLE_ERROR(BF("Could not get velocities from openmm::state"));
	}
	units::Quantity_sp res = units::Quantity_O::create(vel_in_nm_per_ps,
							   _lisp->symbol(units::_sym_nm_per_ps)->dyn<units::Unit_O>(),
							   _lisp);
        return res;
    };

#if 0
    core::CoordinateArray_sp State_O::getForces() const
    {_OF();
	// Get the forces in whatever units they come and convert them to SI 
	core::CoordinateArray_sp res;
	try
	{
	    res = translate::to_object<const std::vector< OpenMM::Vec3 >& >::convert(_lisp,this->_State.getForces());
	}
	catch (...)
	{
	    SIMPLE_ERROR(BF("Could not get forces from openmm::state"));
	}
        return res;
    };
#endif




    units::Quantity_sp State_O::getKineticEnergy()
    {_OF();
	double e_in_kj_per_mole;
	try
	{
	    e_in_kj_per_mole = this->_State.getKineticEnergy();
	}
	catch (...)
	{
	    SIMPLE_ERROR(BF("Could not get kinetic energy from openmm::state"));
	}
	units::Quantity_sp f =
	    units::Quantity_O::create(e_in_kj_per_mole,
				      _lisp->symbol(units::_sym_kilojoules_per_mole)->dyn<units::Unit_O>(),
				      _lisp);
	return f;
    }

    units::Quantity_sp State_O::getPotentialEnergy()
    {_OF();
	double e_in_kj_per_mole;
	try
	{
	     e_in_kj_per_mole = this->_State.getPotentialEnergy();
	}
	catch (...)
	{
	    SIMPLE_ERROR(BF("Could not get potential energy from openmm::state"));
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
{_OF();
    OpenMM::Vec3 a,b,c;
    this->_State.getPeriodicBoxVectors(a,b,c);
    core::Cons_sp res = core::Cons_O::createList(translate::to_object<OpenMM::Vec3>::convert(_lisp,a),
						 translate::to_object<OpenMM::Vec3>::convert(_lisp,b),
						 translate::to_object<OpenMM::Vec3>::convert(_lisp,c),
				   _lisp);
    return res;
}
#endif


    EXPOSE_CLASS(omm,State_O);


};
