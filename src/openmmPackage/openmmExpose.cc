/*
    File: openmmExpose.cc
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
#define DEBUG_LEVEL_NONE


#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/cons.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/package.h>
#include <clasp/core/stringList.h>
#include <core/renderController.h>
#include <clasp/core/lambdaListHandler.h>
#include <clasp/core/numbers.h>
#include <cando/chem/atom.h>
#include <cando/units/unitsPackage.h>
#include <cando/units/quantity.h>
#include <cando/chem/bond.h>
#include <cando/chem/angle.h>
#include <cando/chem/properTorsion.h>
#include <cando/chem/improperTorsion.h>
#include <cando/openmmPackage/openmmPackage.h>
#include <cando/openmmPackage/openmmState.h>
#include <cando/openmmPackage/openmmExpose.h>
#include <external_wrappers.h>
#include <cando/chem/forceField.h>
#include <cando/chem/ffPtorDb.h>
#include <wrappers.h>

namespace omm
{


#define omm_StorePredefinedSymbols
#include <openmmExpose.scrape.inc>

#define omm_code
#include <openmmExpose.scrape.inc>
#undef omm_code
};
namespace omm
{
#define omm_expose
#include <openmmExpose.scrape.inc>
#undef omm_expose
};



namespace omm
{


#define omm_globalInitCode
#include <openmmExpose.scrape.inc>



#define	GLOBAL_INT(sym,val)						\
    {									\
	core::Symbol_sp ss = lisp->internWithPackageName(OmmPkg,#sym)->exportYourself(); \
	ss->defconstant(core::DoubleFloat_O::create(val)); \
    }

    void initializeOpenmmConstants(core::LispPtr lisp)
    {
	LOG("Initializing constants and enums for OpenMM");
	lisp->inPackage(OmmPkg);

#define omm_CreatePredefinedSymbols
#include <openmmExpose.scrape.inc>

#define omm_constants
#include <openmmExpose.scrape.inc>
#undef omm_constants

#define omm_Functions
#include <openmmExpose.scrape.inc>

    }



    template <>
    core::Instance_sp metaClassForOmmClass<OpenMM::Platform>(core::LispPtr lisp)
    {
	return lisp->classFromClassSID(Platform_O::static_classSID());
    }

    template <>
    core::Instance_sp metaClassForOmmClass<OpenMM::LangevinIntegrator>(core::LispPtr lisp)
    {
	return lisp->classFromClassSID(LangevinIntegrator_O::static_classSID());
    }

    template <>
    core::Instance_sp metaClassForOmmClass<OpenMM::VerletIntegrator>(core::LispPtr lisp)
    {
	return lisp->classFromClassSID(VerletIntegrator_O::static_classSID());
    }


    template <>
    core::Instance_sp metaClassForOmmClass<OpenMM::System>(core::LispPtr lisp)
    {
	return lisp->classFromClassSID(System_O::static_classSID());
    }


    
    core::ExternalObject_sp wrap_ommObject(void* ptr, core::Instance_sp mc, core::LispPtr lisp)
    {
	if ( ptr == NULL ) return lisp->nil<core::ExternalObject_O>();
	if ( !lisp->hiddenBinder()->contains(_sym_OmmObjectWrappers)) 
	{
	    lisp->hiddenBinder()->extend(_sym_OmmObjectWrappers,
					 core::ExternalObjectManager_O::create(lisp));
	}
	core::ExternalObjectManager_sp manager =
	    lisp->hiddenBinder()->lookup(_sym_OmmObjectWrappers).as<core::ExternalObjectManager_O>();
	if ( manager->recognizesExternal(ptr) )
	{
	    core::ExternalObject_sp oobj = manager->objectForExternal(ptr);
	    if ( oobj->notNil() ) return oobj;
	}
	core::ExternalObject_sp wrapper = mc->allocate_newNil().as<core::ExternalObject_O>();
	wrapper->initialize();
	wrapper->set_externalObject(ptr);
	manager->registerExternal(ptr,wrapper,lisp);
	return wrapper;
    }


    void System_O::takeOwnershipOfSystem(OpenMM::System* sys)
    {
	if ( this->_externalSystem) delete this->_externalSystem;
	this->_externalSystem = sys;
    }



    core::Cons_sp System_O::getDefaultPeriodicBoxVectors()
    {
	OpenMM::Vec3 a, b, c;
	this->wrappedPtr()->getDefaultPeriodicBoxVectors(a,b,c);
	core::Cons_sp result = core::Cons_O::createList(translate::to_object<OpenMM::Vec3>::convert(_lisp,a),
							translate::to_object<OpenMM::Vec3>::convert(_lisp,b),
							translate::to_object<OpenMM::Vec3>::convert(_lisp,c),
					   this->lisp());
	return result;
    }

    Force_sp System_O::getForce(int index)
    {
	if ( index <0 ) 
	{
	    SIMPLE_ERROR("Index must be positive");
	}
	if ( index >= (int)this->wrappedPtr()->getNumForces() )
	{
	    SIMPLE_ERROR("Illegal index[{}] for force - must be less than[{}]"
			       , index , this->wrappedPtr()->getNumForces() );
	}
	OpenMM::Force* rawForce = &(this->wrappedPtr()->getForce(index));
        core::clasp_write_string(fmt::format("Getting force[{}] = {}\n" , index , rawForce));
	return Force_O::nil(_lisp);
    }


    units::Quantity_sp System_O::getParticleMass(int i)
    {
	double amu = this->wrappedPtr()->getParticleMass(i);
	return units::Quantity_O::create(amu,_lisp->symbol(units::_sym_daltons)->dyn<units::Unit_O>(),_lisp);
    }



    int System_O::addRawForce(OpenMM::Force* rawForce)
    {
	int idx = this->wrappedPtr()->addForce(rawForce);
        core::clasp_write_string(fmt::format("Added force[{}] = {}\n" , idx , rawForce ));
	return idx;
    }


    int System_O::addForce(Force_sp force)
    {
	OpenMM::Force* rawForce = static_cast<OpenMM::Force*>(force->externalObject());
	return this->addRawForce(rawForce);
    }


    /*! Add an Atom to the System and set the atoms ParticleIndex property to the index */
    int	System_O::addCandoAtom(chem::Atom_sp atom)
    {
	double mass = atom->getAtomicWeight();
	int index = this->wrappedPtr()->addParticle(mass);
	atom->setProperty(_lisp->internWithPackageName(Pkg(),ParticleIndex),core::clasp_make_fixnum(index));
	return index;
    }

    core::Cons_sp System_O::getConstraintParameters(int index)
    {
	int particle1, particle2;
	double distance;
	this->wrappedPtr()->getConstraintParameters(index,particle1,particle2,distance);
	core::Cons_sp res = core::Cons_O::createList(core::clasp_make_fixnum(particle1),
					core::clasp_make_fixnum(particle2),
					core::DoubleFloat_O::create(distance),
					_lisp);
	return res;
    }



    Context_sp Context_O::create(System_sp system, Integrator_sp integrator, core::LispPtr lisp)
    {
	Context_sp context = RP_Create<Context_O>(lisp);
	context->setup(system,integrator,lisp->nil<Platform_O>());
	return context;
    }

    void Context_O::setup(System_sp system, Integrator_sp integrator, Platform_sp platform )
    {
	if ( platform.nilp() )
	{
	    this->_externalContext = new OpenMM::Context(*system->wrappedPtr(),
						 *(static_cast<OpenMM::Integrator*>(integrator->externalObject())));
	} else
	{
	    this->_externalContext = new OpenMM::Context(*system->wrappedPtr(),
						 *(static_cast<OpenMM::Integrator*>(integrator->externalObject())),
						 *platform->wrappedPtr() );
	}
    }


    Platform_sp Context_O::getPlatform()
    {
	core::Instance_sp mc = metaClassForOmmClass<OpenMM::Platform>(_lisp);
	Platform_sp platform = wrap_ommObject(&(this->_externalContext->getPlatform()),mc,_lisp).as<Platform_O>();
	return platform;
    }

    System_sp Context_O::getSystem() const
    {
	core::Instance_sp mc = metaClassForOmmClass<OpenMM::System>(_lisp);
	System_sp system = wrap_ommObject(&(this->_externalContext->getSystem()),mc,_lisp).as<System_O>();
	return system;
    }
	

    Integrator_sp Context_O::getIntegrator() const
    {
	OpenMM::Integrator* integrator = &(this->_externalContext->getIntegrator());
	{
	    OpenMM::LangevinIntegrator* subInt = dynamic_cast<OpenMM::LangevinIntegrator*>(integrator);
	    if ( subInt != NULL )
	    {	
		core::Instance_sp mc = metaClassForOmmClass<OpenMM::LangevinIntegrator>(_lisp);
		LangevinIntegrator_sp result = wrap_ommObject(subInt,mc,_lisp).as<LangevinIntegrator_O>();
		return result;
	    }
	}
	{
	    OpenMM::VerletIntegrator* subInt = dynamic_cast<OpenMM::VerletIntegrator*>(integrator);
	    if ( subInt != NULL )
	    {	
		core::Instance_sp mc = metaClassForOmmClass<OpenMM::VerletIntegrator>(_lisp);
		VerletIntegrator_sp result = wrap_ommObject(subInt,mc,_lisp).as<VerletIntegrator_O>();
		return result;
	    }
	}
	SIMPLE_ERROR("Unknown integrator");
    }
	

    ARGUMENT_SC_(positions);
    ARGUMENT_SC_(velocities);
    ARGUMENT_SC_(forces);
    ARGUMENT_SC_(energy);
    ARGUMENT_SC_(parameters);
    ARGUMENT_SC_(enforcePeriodicBox);
    core::T_sp Context_O::getState(core::Function_sp exec, core::Cons_sp args,
				  core::Environment_sp env, core::LispPtr lisp)
    {
	State_sp state = State_O::create(_lisp);
	int types;
	bool positions = env->lookup(_lisp->symbol(_sym_positions))->isTrue();
	bool velocities = env->lookup(_lisp->symbol(_sym_velocities))->isTrue();
	bool forces = env->lookup(_lisp->symbol(_sym_forces))->isTrue();
	bool energy = env->lookup(_lisp->symbol(_sym_energy))->isTrue();
	bool parameters = env->lookup(_lisp->symbol(_sym_parameters))->isTrue();
	bool enforcePeriodicBox = env->lookup(_lisp->symbol(_sym_enforcePeriodicBox))->isTrue();
	types = 0;
	types |= positions ? OpenMM::State::Positions : 0;
	types |= velocities ? OpenMM::State::Velocities : 0;
	types |= forces ? OpenMM::State::Forces : 0;
	types |= energy ? OpenMM::State::Energy : 0;
	types |= parameters ? OpenMM::State::Parameters : 0;
	state->_State = this->wrappedPtr()->getState(types,enforcePeriodicBox);
	return state;
    }


    void Context_O::takeOwnershipOfContext(OpenMM::Context* context)
    {
	if ( this->_externalContext ) delete this->_externalContext;
	this->_externalContext = context;
    }


    void Context_O::setPositions(units::Quantity_sp coords)
    {
	vector<OpenMM::Vec3> positions;
	positions.resize(coords->size());
	for ( int i=0; i<(int)coords->size(); i++ )
	{
	    Vector3 orig = coords->getElement_in_unit_asVector3(i,_lisp->symbol(units::_sym_nanometers)->dyn<units::Unit_O>());
	    OpenMM::Vec3 vec(orig.getX(),orig.getY(),orig.getZ());
	    positions[i] = vec;
	}
	this->wrappedPtr()->setPositions(positions);
    }

    void Context_O::setVelocities(units::Quantity_sp coords)
    {
	ASSERTF(coords->notNil(),("velocities coordinate array is nil"));
	ASSERTF(coords->size()!=0,("There are no velocities"));
	vector<OpenMM::Vec3> velocities;
	velocities.resize(coords->size());
	for ( int i=0; i<(int)coords->size(); i++ )
	{
	    Vector3 orig = coords->getElement_in_unit_asVector3(i,_lisp->symbol(units::_sym_nm_per_ps)->dyn<units::Unit_O>());
	    OpenMM::Vec3 vec(orig.getX(),orig.getY(),orig.getZ());
	    velocities[i] = vec;
	}
	this->wrappedPtr()->setVelocities(velocities);
    }


    int HarmonicBondForce_O::addCandoBond(chem::Bond_sp b, chem::ForceField_sp forceField)
    {
	chem::Atom_sp a1 = b->getAtom1();
	chem::Atom_sp a2 = b->getAtom2();
	chem::FFStretch_sp term = forceField->getStretchDb()->findTerm(a1,a2);
	double r0 = term->getR0_Nanometer();
	double kb = term->getKb_kjPerNanometerSquared();
	ASSERT(_Force!=NULL);
	int particle1 = a1->getProperty(_lisp->internWithPackageName(Pkg(),ParticleIndex)).as<core::Fixnum_O>()->get();
	int particle2 = a2->getProperty(_lisp->internWithPackageName(Pkg(),ParticleIndex)).as<core::Fixnum_O>()->get();
	int bondId = _Force->addBond(particle1,particle2,r0,kb);
	return bondId;
    };



    core::Cons_sp HarmonicBondForce_O::getBondParameters(int index) const
    {
	int particle1, particle2;
	double length, k;
	ASSERT(this->_Force!=NULL);
	this->_Force->getBondParameters(index,particle1,particle2,length,k);
	core::Cons_sp result = core::Cons_O::createList(
	    core::clasp_make_fixnum(particle1),
	    core::clasp_make_fixnum(particle2),
	    core::DoubleFloat_O::create(length),
	    core::DoubleFloat_O::create(k),_lisp);
	return result;
    };



    int HarmonicAngleForce_O::addCandoAngle(chem::Angle_sp ang, chem::ForceField_sp forceField)
    {
	chem::Atom_sp a1 = ang->getA1();
	chem::Atom_sp a2 = ang->getA2();
	chem::Atom_sp a3 = ang->getA3();
	chem::FFAngle_sp term = forceField->getAngleDb()->findTerm(a1,a2,a3);
	double t0 = term->getAngle_Radian();
	double kt = term->getK2_kjPerRadianSquared();
	ASSERT(_Force!=NULL);
	int particle1 = a1->getProperty(_lisp->internWithPackageName(Pkg(),ParticleIndex)).as<core::Fixnum_O>()->get();
	int particle2 = a2->getProperty(_lisp->internWithPackageName(Pkg(),ParticleIndex)).as<core::Fixnum_O>()->get();
	int particle3 = a3->getProperty(_lisp->internWithPackageName(Pkg(),ParticleIndex)).as<core::Fixnum_O>()->get();
	int angleId = _Force->addAngle(particle1,particle2,particle3,t0,kt);
	return angleId;
    };



    core::Cons_sp HarmonicAngleForce_O::getAngleParameters(int index) const
    {
	int particle1, particle2, particle3;
	double angle, k;
	ASSERT(this->_Force!=NULL);
	this->_Force->getAngleParameters(index,particle1,particle2,particle3,angle,k);
	core::Cons_sp result = core::Cons_O::createList(
	    core::clasp_make_fixnum(particle1),
	    core::clasp_make_fixnum(particle2),
	    core::clasp_make_fixnum(particle3),
	    core::DoubleFloat_O::create(angle),
	    core::DoubleFloat_O::create(k),_lisp);
	return result;
    };



    //
    // PeriodicTorsionForce
    //


    /*! Add all of the Torsion terms for each periodicity in the forceField
      return a Cons of indexes for the new Torsion terms added */
    core::Cons_sp PeriodicTorsionForce_O::addCandoProperTorsions(chem::ProperTorsion_sp tor,
								 chem::ForceField_sp forceField)
    {
	chem::Atom_sp a1 = tor->getA1();
	chem::Atom_sp a2 = tor->getA2();
	chem::Atom_sp a3 = tor->getA3();
	chem::Atom_sp a4 = tor->getA4();
	string t1 = a1->getTypeString();
	string t2 = a2->getTypeString();
	string t3 = a3->getTypeString();
	string t4 = a4->getTypeString();
	if ( !forceField->getPtorDb()->hasBestTerm(t1,t2,t3,t4) )
	{
	    SIMPLE_ERROR("Could not find torsion term for [{}-{}-{}-{}]"
			       , t1 , t2 %t3 %t4 );
	}
	chem::FFPtor_sp ffPtor = forceField->getPtorDb()->findBestTerm(t1,t2,t3,t4);
	core::Symbol_sp particleIndexSymbol = _lisp->internWithPackageName(Pkg(),ParticleIndex);
	core::Cons_sp first = _lisp->create<core::Cons_O>(core::Cons_O::nil(_lisp));
	core::Cons_sp cur = first;
	for ( int period=1; period<=chem::FFPtor_O::MaxPeriodicity; period++ )
	{
	    if ( ffPtor->hasPeriodicity(period) )
	    {
		LOG("Adding Proper Torsion period[{}] {}/{}-{}/{}-{}/{}-{}/{}"
		    , period
		    , a1->getName() , t1
		    , a2->getName() , t2
		    , a3->getName() , t3
		    , a4->getName() , t4 );
		int i1 = a1->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
		int i2 = a2->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
		int i3 = a3->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
		int i4 = a4->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
		int idx = this->wrappedPtr()->addTorsion(i1,i2,i3,i4,period,ffPtor->getPhaseRadians(period),
					    ffPtor->getV_kj(period));
		core::Cons_sp one = _lisp->create<core::Cons_O>(core::clasp_make_fixnum(idx));
		cur->setCdr(one);
		cur = one;
	    }
	}
	return first->cdr();
    }

    core::Cons_sp PeriodicTorsionForce_O::addCandoImproperTorsions(chem::ImproperTorsion_sp tor,
								   chem::ForceField_sp forceField)
    {
	chem::Atom_sp a1 = tor->getA1();
	chem::Atom_sp a2 = tor->getA2();
	chem::Atom_sp a3central = tor->getA3Central();
	chem::Atom_sp a4 = tor->getA4();
	string t1 = a1->getTypeString();
	string t2 = a2->getTypeString();
	string t3 = a3central->getTypeString();
	string t4 = a4->getTypeString();
	core::Symbol_sp particleIndexSymbol = _lisp->internWithPackageName(Pkg(),ParticleIndex);
	core::Cons_sp first = _lisp->create<core::Cons_O>(core::Cons_O::nil(_lisp));
	core::Cons_sp cur = first;
	if ( forceField->getItorDb()->hasBestTerm(t1,t2,t3,t4) ) 
	{
	    chem::FFItor_sp ffItor = forceField->_Itors->findBestTerm(t1,t2,t3,t4);
	    for ( int n=1;n<=chem::FFItor_O::IMaxPeriodicity; n++ )
	    {
		if ( ffItor->hasPeriodicity(n) )
		{
		    LOG("Adding Improper Torsion period[{}] {}/{}-{}/{}-{}/{}-{}/{}"
			% n
			% a1->getName() , t1
			% a2->getName() , t2
			% a3central->getName() , t3
			% a4->getName() , t4 );
		    int i1 = a1->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
		    int i2 = a2->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
		    int i3 = a3central->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
		    int i4 = a4->getProperty(particleIndexSymbol).as<core::Fixnum_O>()->get();
		    int idx = this->wrappedPtr()->addTorsion(i1,i2,i3,i4,n,ffItor->getPhaseRadians(n),
							  ffItor->getV_kj(n));
		    core::Cons_sp one = _lisp->create<core::Cons_O>(core::clasp_make_fixnum(idx));
		    cur->setCdr(one);
		    cur = one;
		}
	    }
	}
	return first->cdr();
    }


    core::Cons_sp PeriodicTorsionForce_O::getTorsionParameters(int idx) const
    {
	IMPLEMENT_ME();
    }

    // 
    // NonbondedForce

    void NonbondedForce_O::initialize()
    {
	this->Base::initialize();
    }


    int NonbondedForce_O::addNonbondedAtom(chem::Atom_sp atom, chem::ForceField_sp forceField)
    {
	string ty = atom->getTypeString();
	if ( !forceField->getNonbondDb()->hasType(ty) )
	{
	    SIMPLE_ERROR("Could not find nonbond type[{}]" , ty );
	}
	double charge = atom->getCharge();
	chem::FFNonbond_sp nonbond = forceField->getNonbondDb()->findType(ty);
	double sigma = nonbond->getRadius_Nanometers();
	double epsilon = nonbond->getEpsilon_kj();
	int idx = this->wrappedPtr()->addParticle(charge,sigma,epsilon);
	int recordedIdx = atom->getProperty(_lisp->internWithPackageName(Pkg(),ParticleIndex)).as<core::Fixnum_O>()->get();
	ASSERTF(idx==recordedIdx,("There was a problem, you added a particle to the NonbondedForce and it returned index[{}] but it should have been the index[{}] for atom[{}] - make sure you add atoms to the NonbondForce immediately after you add them to the System") , idx , recordedIdx , atom->description() );
	return idx;
    }

    core::Symbol_sp NonbondedForce_O::getNonbondedMethod()
    {
	OpenMM::NonbondedForce::NonbondedMethod nbm = this->wrappedPtr()->getNonbondedMethod();
	switch (nbm)
	{
	case OpenMM::NonbondedForce::NoCutoff:
	{
	    return _sym_NoCutoff;
	}
	case OpenMM::NonbondedForce::CutoffNonPeriodic:
	{
	    return _sym_CutoffNonPeriodic;
	}
	case OpenMM::NonbondedForce::CutoffPeriodic:
	{
	    return _sym_CutoffPeriodic;
	}
	case OpenMM::NonbondedForce::Ewald:
	{
	    return _sym_Ewald;
	}
	case OpenMM::NonbondedForce::PME:
	{
	    return _sym_PME;
	}
	}
	SIMPLE_ERROR("Illegal NonbondedMethod");
    }








    void NonbondedForce_O::setNonbondedMethod(core::Symbol_sp sym)
    {
	OpenMM::NonbondedForce::NonbondedMethod nbm;
	if ( sym == _sym_NoCutoff )
	{
	    nbm = OpenMM::NonbondedForce::NoCutoff;
	} else if ( sym == _sym_CutoffNonPeriodic )
	{
	    nbm = OpenMM::NonbondedForce::CutoffNonPeriodic;
	} else if ( sym == _sym_CutoffPeriodic )
	{
	    nbm = OpenMM::NonbondedForce::CutoffPeriodic;
	} else if ( sym == _sym_Ewald )
	{
	    nbm = OpenMM::NonbondedForce::Ewald;
	} else if ( sym == _sym_PME)
	{
	    nbm = OpenMM::NonbondedForce::PME;
	} else
	{
          SIMPLE_ERROR("Illegal NonbondedMethod[{}]" , _rep_(sym) );
	}
	this->wrappedPtr()->setNonbondedMethod(nbm);
    }

    core::Cons_sp NonbondedForce_O::getParticleParameters(int p)
    {
	double charge, sigma, epsilon;
	this->wrappedPtr()->getParticleParameters(p,charge,sigma,epsilon);
	core::Cons_sp res = core::Cons_O::createList(core::DoubleFloat_O::create(charge),
						     core::DoubleFloat_O::create(sigma),
						     core::DoubleFloat_O::create(epsilon),_lisp);
	return res;
    }

    core::Cons_sp NonbondedForce_O::getExceptionParameters(int p)
    {
	int part1, part2;
	double chargeProd, sigma, epsilon;
	this->wrappedPtr()->getExceptionParameters(p,part1,part2,chargeProd,sigma,epsilon);
	core::Cons_sp res = core::Cons_O::createList(core::clasp_make_fixnum(part1),
						     core::clasp_make_fixnum(part2),
						     core::DoubleFloat_O::create(chargeProd),
						     core::DoubleFloat_O::create(sigma),
						     core::DoubleFloat_O::create(epsilon),_lisp);
	return res;
    }

    int NonbondedForce_O::addException(int part1, int part2, double chargeProd, double sigma, double epsilon, bool replace )
    {
	int idx = this->wrappedPtr()->addException(part1,part2,chargeProd,sigma,epsilon,replace);
	return idx;
    }



    void VerletIntegrator_O::initialize()
    {
	this->Base::initialize();
    }




    int wrap_Platform_getNumPlatforms(core::LispPtr lisp)
    {
	return OpenMM::Platform::getNumPlatforms();
    }

    Platform_sp wrap_Platform_getPlatform(int idx, core::LispPtr lisp)
    {
	if ( idx < 0 || idx >= OpenMM::Platform::getNumPlatforms() ) 
	{
	    THROW(lisp->error(BF("Illegal platform index[%d] must be between 0 and %d") , idx , (OpenMM::Platform::getNumPlatforms()-1) ));
	}
	OpenMM::Platform& p = OpenMM::Platform::getPlatform(idx);
	Platform_sp op = wrap_ommObject(&p,metaClassForOmmClass<OpenMM::Platform>(lisp),lisp).as<Platform_O>();
	return op;
    }

    Platform_sp wrap_Platform_getPlatformByName(string const& name, core::LispPtr lisp)
    {
	OpenMM::Platform& p = OpenMM::Platform::getPlatformByName(name);
	Platform_sp op = wrap_ommObject(&p,metaClassForOmmClass<OpenMM::Platform>(lisp),lisp).as<Platform_O>();
	return op;
    }


    core::Cons_sp wrap_Platform_loadPluginsFromDirectory( string const& directoryName, core::LispPtr lisp)
    {
	vector<string> plugins = OpenMM::Platform::loadPluginsFromDirectory(directoryName);
	core::Cons_sp first = core::Cons_O::create(_lisp->onil(),_lisp);
	core::Cons_sp cur = first;
	for ( vector<string>::iterator si=plugins.begin(); si!=plugins.end(); si++ )
	{
	    core::Cons_sp one = _lisp->create<core::Cons_O>(core::Str_O::create(*si,_lisp));
	    cur->setCdr(one);
	    cur = one;
	}
	return first->cdr();
    }


};
