include(`exposeMacros.m4')

#ifdef	OMM_NAMESPACE_DECLARATION
namespace omm
{
#endif

// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------
//
// Bindings
//
// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------


ommGlobalInt(Positions,OpenMM::State::Positions);
ommGlobalInt(Velocities,OpenMM::State::Velocities);
ommGlobalInt(Forces,OpenMM::State::Forces);
ommGlobalInt(Parameters,OpenMM::State::Parameters);
ommGlobalInt(Energy,OpenMM::State::Energy);

ommSymbol(OmmObjectWrappers)




// By default we expand into the omm_header block
//
#ifdef omm_header


SYMBOL_SC_(Context);
ommClass(Context,core::ExternalObject_O)
 wrappedObject(OpenMM::Context,_externalContext);
     Context_O(const core::MetaClass_sp& mc) : Base(mc), _externalContext(NULL) {};
     virtual ~Context_O() { if ( this->_externalContext) delete this->_externalContext; };
    lispInitialization("(self &key omm::system omm::integrator omm::platform)",`
        System_sp system = bargs->lookup(OmmPkg,"system")->as<System_O>();
	Integrator_sp integrator = bargs->lookup(OmmPkg,"integrator")->as<Integrator_O>();
	Platform_sp platform = bargs->lookup(OmmPkg,"platform")->as<Platform_O>();
	this->setup(system,integrator,platform);
	return lisp->onil();
    ')
public:
	static Context_sp create(System_sp system, Integrator_sp integrator, core::Lisp_sp lisp);
public:
	void setup(System_sp system, Integrator_sp integrator, Platform_sp platform );
 oMethodDontExpose(takeOwnershipOfContext,`void takeOwnershipOfContext(OpenMM::Context* context)')
    core::T_sp getState(core::Function_sp exec, core::Cons_sp args,
				  core::Environment_sp env, core::Lisp_sp lisp);
 noWrapMethodWithArguments(getState,(self &key positions velocities forces energy parameters enforce-periodic-box))
 oMethod(getPlatform,`Platform_sp getPlatform();')
 ommMethod(setTime)
 oMethod(setPositions,`void setPositions(units::Quantity_sp positions)')
 oMethod(setVelocities,`void setVelocities(units::Quantity_sp velocities)')
 oMethod(getSystem,`System_sp getSystem() const;')
 oMethod(getIntegrator,`Integrator_sp getIntegrator() const;')
 ommMethod(getParameter)
 ommMethod(setParameter)
 ommMethod(setPeriodicBoxVectors)
 ommMethod(applyConstraints)
 ommMethod(reinitialize)
end_ommClass_noDefaultCtorDtor()





SYMBOL_SC_(Force);
ommClass(Force,core::ExternalObject_O)
end_ommClass()

SYMBOL_SC_(HarmonicBondForce);
ommClass(HarmonicBondForce,Force_O)
 wrappedObject(OpenMM::HarmonicBondForce,_Force)
 public:
     HarmonicBondForce_O(const core::MetaClass_sp& mc) : Base(mc), _Force(NULL)
     { _Force = new OpenMM::HarmonicBondForce();};
    lispInitialization("(self omm::system)",`
        System_sp system = bargs->lookup(Pkg(),"system")->as<System_O>();
	system->addRawForce(this->_Force);
	return lisp->onil();
    ')
 ommMethod(getNumBonds)
 ommMethod(addBond)
 oMethod(addCandoBond,`int addCandoBond(chem::Bond_sp b, chem::ForceField_sp forceField)')
 ommMethod(setBondParameters)
 oMethod(getBondParameters,`core::Cons_sp getBondParameters(int index) const')
end_ommClass_noDefaultCtorDtor()



SYMBOL_SC_(HarmonicAngleForce);
ommClass(HarmonicAngleForce,Force_O)
 wrappedObject(OpenMM::HarmonicAngleForce,_Force)
 public:
     HarmonicAngleForce_O(const core::MetaClass_sp& mc) : Base(mc), _Force(NULL)
          { _Force = new OpenMM::HarmonicAngleForce();};
    lispInitialization("(self omm::system)",`
        System_sp system = bargs->lookup(Pkg(),"system")->as<System_O>();
	system->addRawForce(this->_Force);
	return lisp->onil();
    ')
 ommMethod(getNumAngles)
 ommMethod(addAngle)
 oMethod(addCandoAngle,`int addCandoAngle(chem::Angle_sp ang, chem::ForceField_sp forceField)')
 ommMethod(setAngleParameters)
 oMethod(getAngleParameters,`core::Cons_sp getAngleParameters(int index) const')
end_ommClass_noDefaultCtorDtor()



SYMBOL_SC_(PeriodicTorsionForce);
ommClass(PeriodicTorsionForce,Force_O)
 wrappedObject(OpenMM::PeriodicTorsionForce,_Force)
 public:
     PeriodicTorsionForce_O(const core::MetaClass_sp& mc) : Base(mc), _Force(NULL)
          { _Force = new OpenMM::PeriodicTorsionForce();};
    lispInitialization("(self omm::system)",`
        System_sp system = bargs->lookup(Pkg(),"system")->as<System_O>();
	system->addRawForce(this->_Force);
	return lisp->onil();
    ')
 ommMethod(getNumTorsions)
 ommMethod(addTorsion)
 oMethod(addCandoProperTorsions,`core::Cons_sp addCandoProperTorsions(chem::ProperTorsion_sp tor, chem::ForceField_sp forceField)')
 oMethod(addCandoImproperTorsions,`core::Cons_sp addCandoImproperTorsions(chem::ImproperTorsion_sp tor, chem::ForceField_sp forceField)')
 ommMethod(setTorsionParameters)
 oMethod(getTorsionParameters,`core::Cons_sp getTorsionParameters(int index) const')
end_ommClass_noDefaultCtorDtor()








SYMBOL_SC_(NonbondedForce);
ommClass(NonbondedForce,Force_O)
 oMethodDontExpose(initialize,`void initialize();')
 wrappedObject(OpenMM::NonbondedForce,_Force)
 public:
     NonbondedForce_O(const core::MetaClass_sp& mc) : Base(mc), _Force(NULL)
     { _Force = new OpenMM::NonbondedForce();};
 lispInitialization("(self omm::system)",`
        System_sp system = bargs->lookup(Pkg(),"system")->as<System_O>();
	system->addRawForce(this->_Force);
	return lisp->onil();
    ')
 ommMethod(getNumParticles)
 ommMethod(getNumExceptions)
 ommSymbol(NoCutoff)
 ommSymbol(CutoffNonPeriodic)
 ommSymbol(CutoffPeriodic)
 ommSymbol(Ewald)
 ommSymbol(PME)
 oMethod(getNonbondedMethod,`core::Symbol_sp getNonbondedMethod();')
 oMethod(setNonbondedMethod,`void setNonbondedMethod(core::Symbol_sp sym);')
 ommMethod(getCutoffDistance)
 ommMethod(setCutoffDistance)
 ommMethod(getReactionFieldDielectric)
 ommMethod(setReactionFieldDielectric)
 ommMethod(getEwaldErrorTolerance)
 ommMethod(setEwaldErrorTolerance)
 ommMethod(addParticle)
 oMethod(addNonbondedAtom,`int addNonbondedAtom(chem::Atom_sp a, chem::ForceField_sp forceField)')
 oMethod(getParticleParameters,`core::Cons_sp getParticleParameters(int p);')
 ommMethod(setParticleParameters)
 oMethod(addException,`int addException(int particle1, int particle2, double chargeProd, double sigma, double epsilon, bool replace);')
 oMethod(getExceptionParameters,`core::Cons_sp getExceptionParameters(int p);')
 ommMethod(setExceptionParameters)
 ommMethod(getUseDispersionCorrection)
 ommMethod(setUseDispersionCorrection)
end_ommClass_noDefaultCtorDtor()


SYMBOL_SC_(Integrator);
ommClass(Integrator,core::ExternalObject_O)
 ommMethod(getStepSize)
 ommMethod(setStepSize)
 ommMethod(getConstraintTolerance)
 ommMethod(setConstraintTolerance)
 ommMethod(step)
end_ommClass()

SYMBOL_SC_(LangevinIntegrator);
ommClass(LangevinIntegrator,Integrator_O)
 wrappedObject(OpenMM::LangevinIntegrator,_Integrator)
 public:
     LangevinIntegrator_O(const core::MetaClass_sp& mc) : Base(mc), _Integrator(NULL) {}
     lispInitialization(`"(self &key (temperature (core:* units:kelvin 300.0))"
     			       	    "(frictionCoeff (core:* (units:Unit units:per_picosecond) 90.0))"
				    "(timeStep (core:* units:femtoseconds 1.0)))"',`
        units::Quantity_sp temperature = bargs->lookup(OmmPkg,"temperature")->as<units::Quantity_O>();
        units::Quantity_sp frictionCoeff = bargs->lookup(OmmPkg,"frictionCoeff")->as<units::Quantity_O>();
        units::Quantity_sp timeStep = bargs->lookup(OmmPkg,"timeStep")->as<units::Quantity_O>();
	{
	    double tempK = temperature->value_in_unit_asReal(_lisp->symbol(units::_sym_kelvin)->dyn<units::Unit_O>());
	    double frictionCoeff_per_ps = frictionCoeff->value_in_unit_asReal(_lisp->symbol(units::_sym_per_picosecond)->dyn<units::Unit_O>());
	    double timeStep_ps = timeStep->value_in_unit_asReal(_lisp->symbol(units::_sym_picoseconds)->dyn<units::Unit_O>());
	    LOG(BF("Creating LangevinIntegrator with temp[%lf] frictionCoeff[%lf] timeStep[%lf]\n") % tempK % frictionCoeff_per_ps %timeStep_ps);
	    this->_Integrator = new OpenMM::LangevinIntegrator(tempK,frictionCoeff_per_ps,timeStep_ps);
	};
	return lisp->onil();
    ')
 ommMethod(getTemperature)
 ommMethod(setTemperature)
 ommMethod(getFriction)
 ommMethod(setFriction)
 ommMethod(getRandomNumberSeed)
 ommMethod(setRandomNumberSeed)
end_ommClass_noDefaultCtorDtor()



SYMBOL_SC_(VerletIntegrator);
ommClass(VerletIntegrator,Integrator_O)
 wrappedObject(OpenMM::VerletIntegrator,_Integrator)
 public:
     VerletIntegrator_O(const core::MetaClass_sp& mc) : Base(mc), _Integrator(NULL) {};
    lispInitialization("(self &key (omm::timeStepFs 0.001))",`
        double timeStepFs = bargs->lookup(OmmPkg,"timeStepFs")->as<core::DoubleFloat_O>()->get();
	{ _Integrator = new OpenMM::VerletIntegrator(timeStepFs);};
	return lisp->onil();
    ')
 oMethodDontExpose(initialize,`void initialize();')
end_ommClass_noDefaultCtorDtor()



SYMBOL_SC_(Platform);
ommClass(Platform,core::ExternalObject_O)
  wrappedObject(OpenMM::Platform,_Platform)
  ommMethod(getName)
  ommMethod(getSpeed)
  ommMethod(supportsDoublePrecision)

  omm_defun(getNumPlatforms,wrap_Platform_getNumPlatforms)
  omm_defun(getPlatform,wrap_Platform_getPlatform)
  omm_defun(getPlatformByName,wrap_Platform_getPlatformByName)
  omm_defun(loadPluginsFromDirectory,wrap_Platform_loadPluginsFromDirectory)
end_ommClass()






SYMBOL_SC_(System);	  
ommClass(System,core::ExternalObject_O)
 wrappedObject(OpenMM::System,_externalSystem)
 public:
     System_O(const core::MetaClass_sp& mc) : Base(mc), _externalSystem(NULL) {this->_externalSystem = new OpenMM::System();};
     virtual ~System_O() { /* if ( this->_externalSystem) delete this->_externalSystem; */ };
 instanceVariable(`core::Vector0<Force_O> _Forces;');
 ommMethod(getNumParticles)
 ommMethod(addParticle)
 oMethodDontExpose(takeOwnershipOfSystem,`void takeOwnershipOfSystem(OpenMM::System* sys)')
 oMethod(addCandoAtom,`int addCandoAtom(chem::Atom_sp atom)')
 oMethod(getParticleMass,`units::Quantity_sp getParticleMass(int i)')
 ommMethod(getNumConstraints)
 ommMethod(addConstraint)
 oMethod(getConstraintParameters,`core::Cons_sp getConstraintParameters(int index)')
 ommMethod(setConstraintParameters)
 oMethod(addForce,`int addForce(Force_sp force);')
 oMethodDontExpose(addRawForce,`int addRawForce(OpenMM::Force* force);')
 ommMethod(getNumForces)
 oMethod(getForce,`Force_sp getForce(int index);')
 oMethod(getDefaultPeriodicBoxVectors,`core::Cons_sp getDefaultPeriodicBoxVectors();')
 ommMethod(setDefaultPeriodicBoxVectors)
end_ommClass_noDefaultCtorDtor()





#endif // close omm_header block


// #ifdef	omm_header
// undivert(fHeader)
// #endif


#ifdef	omm_code
undivert(fCode)
#endif // close omm_code block
#ifdef	omm_expose
undivert(fExpose)
#endif // close omm_expose block
#ifdef	omm_constants
#undef	omm_constants
undivert(fInit)
#endif // close omm_constants block

#ifdef	omm_globalInitCode
#undef	omm_globalInitCode
string ommInitializationCode()
{
    	/* Return a block of Lisp code to initialize omm
	*/
    stringstream sInitCode;
    undivert(fInitCode)
    return sInitCode.str();
}
#endif // close omm_globalInitCode block




#ifdef	omm_DeclarePredefinedSymbols
#undef	omm_DeclarePredefinedSymbols
undivert(fDeclarePredefinedSymbols)
#endif


#ifdef	omm_StorePredefinedSymbols
#undef	omm_StorePredefinedSymbols
undivert(fStorePredefinedSymbols)
#endif


#ifdef	omm_CreatePredefinedSymbols
#undef	omm_CreatePredefinedSymbols
undivert(fCreatePredefinedSymbols)
#endif

#ifdef	omm_Functions
#undef  omm_Functions
undivert(fFunctions)
#endif

#ifdef	omm_Translate
#undef	omm_Translate
undivert(fTranslate)
#endif



#ifdef	OMM_NAMESPACE_DECLARATION
}; // namespace omm
#endif
