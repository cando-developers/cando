#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/rotamer.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(kinematics,Rotamer_O);

    void Rotamer_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Rotamer_O>()
//	    .def_raw("__init__",&Rotamer_O::__init__,"(self &key dihedrals sigmas indices (count 0) (probability 0.0))")
	    .def("rotamer-probability",&Rotamer_O::probability)
	;
    }

    void Rotamer_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,Rotamer,"","",_lisp)
//	.def_raw("__init__",&"(self &key dihedrals sigmas indices (count 0) (probability 0.0))")
	    .def("rotamer-probability",&Rotamer_O::probability)
	;
#endif
    }


#if INIT_TO_FACTORIES

#define ARGS_Rotamer_O_make "(dihedrals sigmas indicies probability count)"
#define DECL_Rotamer_O_make ""
#define DOCS_Rotamer_O_make "make Rotamer"
  Rotamer_sp Rotamer_O::make(core::Cons_sp& dihedrals, core::Cons_sp& sigmas, core::Cons_sp& indices, const double probability, const int count)
    {_G();
	GC_ALLOCATE(Rotamer_O, me );
      me->_Count = count;
      ASSERTF(dihedrals->length()==sigmas->length(),
	      BF("Mismatch between #dihedrals[%d] and #sigmas[%d]")
	      % dihedrals->length() % sigmas->length() );
      me->_Probability = probability;
      for ( ; dihedrals.notnilp(); dihedrals = dihedrals->cdr(),
	      sigmas = sigmas->cdr(),
	      indices = indices->cdr() )
	{
	  me->_Dihedrals.push_back(dihedrals->car<core::DoubleFloat_O>()->get());
	  me->_Sigmas.push_back(sigmas->car<core::DoubleFloat_O>()->get());
	  if ( indices->ocar().notnilp() )
	    {
	      me->_Indices.push_back(indices->car<core::Fixnum_O>()->get());
	    } else
	    {
	      me->_Indices.push_back(0);
	    }
	}
	return me;
    };

#else
    core::T_sp Rotamer_O::__init__(core::Function_sp exec,
				  core::Cons_sp args,
				  core::Environment_sp env,
				  core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
	core::Cons_sp dihedrals = translate::from_object<core::Cons_sp>::convert(env->lookup(lisp->symbol(_sym_dihedrals)));
	core::Cons_sp sigmas = translate::from_object<core::Cons_sp>::convert(env->lookup(lisp->symbol(_sym_sigmas)));
	core::Cons_sp indices = translate::from_object<core::Cons_sp>::convert(env->lookup(lisp->symbol(_sym_indices)));
	double probability = translate::from_object<double>::convert(env->lookup(lisp->symbol(_sym_probability)));
	this->_Count = translate::from_object<int>::convert(env->lookup(lisp->symbol(_sym_count)));
	ASSERTF(dihedrals->length()==sigmas->length(),
		BF("Mismatch between #dihedrals[%d] and #sigmas[%d]")
		% dihedrals->length() % sigmas->length() );
	this->_Probability = probability;
	for ( ; dihedrals.notnilp(); dihedrals = dihedrals->cdr(),
		  sigmas = sigmas->cdr(),
		  indices = indices->cdr() )
	{
	    this->_Dihedrals.push_back(dihedrals->car<core::DoubleFloat_O>()->get());
	    this->_Sigmas.push_back(sigmas->car<core::DoubleFloat_O>()->get());
	    if ( indices->ocar().notnilp() )
	    {
		this->_Indices.push_back(indices->car<core::Fixnum_O>()->get());
	    } else
	    {
		this->_Indices.push_back(0);
	    }
	}
	return _Nil<core::T_O>();
    }
#endif


    string Rotamer_O::__repr__() const
    {
	stringstream ss;
	ss << "(";
	ss << this->className();
	ss << " ";
	ss << ":dihedrals (";
	{
	    for ( int i=0; i<(int)this->_Dihedrals.size(); i++ )
	    {
		ss << this->_Dihedrals[i] << " ";
	    }
	}
	ss << ") :sigmas (";
	{
	    for ( int i=0; i<(int)this->_Sigmas.size(); i++ )
	    {
		ss << this->_Sigmas[i] << " ";
	    }
	}
	ss << ") :probability " << this->_Probability << " )";
	return ss.str();
    }


#ifdef XML_ARCHIVE
    void Rotamer_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->archiveVectorDouble("dih",this->_Dihedrals);
	node->archiveVectorDouble("sig",this->_Sigmas);
	node->archiveVectorInt("ind",this->_Indices);
	node->attribute("c",this->_Count);
	node->attribute("p",this->_Probability);
    }
#endif


    void Rotamer_O::initialize()
    {_OF();
        this->Base::initialize();
    }



    EXPOSE_CLASS(kinematics,RotamerSetBase_O);


    void RotamerSetBase_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<RotamerSetBase_O>()
//	.def_raw("__init__",&RotamerSetBase_O::__init__,"(self)")
	    .def("rotamer-set-add-rotamer",&RotamerSetBase_O::addRotamer)
	    .def("rotamer-set-as-cons",&RotamerSetBase_O::asCons)
	    .def("rotamer-set-size",&RotamerSetBase_O::size)
	    .def("rotamer-set-get",&RotamerSetBase_O::get)
	;
    }

    void RotamerSetBase_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,RotamerSetBase,"","",_lisp)
//	.def_raw("__init__",&"(self)")
	    .def("rotamer-set-add-rotamer",&RotamerSetBase_O::addRotamer)
	    .def("rotamer-set-as-cons",&RotamerSetBase_O::asCons)
	    .def("rotamer-set-size",&RotamerSetBase_O::size)
	    .def("rotamer-set-get",&RotamerSetBase_O::get)
	;
#endif
    }


#if INIT_TO_FACTORIES

#define ARGS_RotamerSetBase_O_make "()"
#define DECL_RotamerSetBase_O_make ""
#define DOCS_RotamerSetBase_O_make "make RotamerSetBase"
  RotamerSetBase_sp RotamerSetBase_O::make()
    {_G();
      IMPLEMENT_ME();
    };

#else
    core::T_sp RotamerSetBase_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void RotamerSetBase_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->archiveVector0("rotamers",this->_Rotamers);
    }
#endif


    void RotamerSetBase_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    Rotamer_sp RotamerSetBase_O::addRotamer(Rotamer_sp rotamer)
    {_OF();
	this->_Rotamers.append(rotamer);
	return rotamer;
    }

    Rotamer_sp RotamerSetBase_O::get(int i) const
    {_OF();
	ASSERTF(i>=0 && i<(int)this->_Rotamers.size(),BF("Illegal index for rotamer[%d]") % i);
	return this->_Rotamers[i];
    }
    

    core::Cons_sp RotamerSetBase_O::asCons() const
    {_OF();
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>());
	core::Cons_sp cur = first;
	for ( int i=0; i<(int)this->_Rotamers.size(); i++ )
	{
	    core::Cons_sp one = core::Cons_O::create(this->_Rotamers[i]);
	    cur->setCdr(one);
	    cur = one;
	}
	return first->cdr();
    }















    EXPOSE_CLASS(kinematics,RotamerSet_O);


    void RotamerSet_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<RotamerSet_O>()
	    ;
    }

    void RotamerSet_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,RotamerSet,"","",_lisp)
	    ;
#endif
    }


#if INIT_TO_FACTORIES

#define ARGS_RotamerSet_O_make "(RotamerSet)"
#define DECL_RotamerSet_O_make ""
#define DOCS_RotamerSet_O_make "make RotamerSet"
  RotamerSet_sp RotamerSet_O::make()
    {_G();
      IMPLEMENT_ME();
    };

#else
    core::T_sp RotamerSet_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void RotamerSet_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
    }
#endif


    void RotamerSet_O::initialize()
    {_OF();
        this->Base::initialize();
    }









    EXPOSE_CLASS(kinematics,BackboneDependentRotamerSet_O);
    void BackboneDependentRotamerSet_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<BackboneDependentRotamerSet_O>()
//	    .def_raw("__init__",&BackboneDependentRotamerSet_O::__init__,"(self &key phi psi)")
	    .def("backbone-dependent-rotamer-set-phi",&BackboneDependentRotamerSet_O::phi)
	    .def("backbone-dependent-rotamer-set-psi",&BackboneDependentRotamerSet_O::psi)
	    ;
    }

    void BackboneDependentRotamerSet_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,BackboneDependentRotamerSet,"","",_lisp)
//	    .def_raw("__init__",&"(self &key phi psi)")
	    .def("backbone-dependent-rotamer-set-phi",&BackboneDependentRotamerSet_O::phi)
	    .def("backbone-dependent-rotamer-set-psi",&BackboneDependentRotamerSet_O::psi)
	    ;
#endif
    }

    void BackboneDependentRotamerSet_O::initialize()
    {
	this->Base::initialize();
    }

#if INIT_TO_FACTORIES

#define ARGS_BackboneDependentRotamerSet_O_make "(phi psi)"
#define DECL_BackboneDependentRotamerSet_O_make ""
#define DOCS_BackboneDependentRotamerSet_O_make "make BackboneDependentRotamerSet"
  BackboneDependentRotamerSet_sp BackboneDependentRotamerSet_O::make(const int phi, const int psi)
    {_G();
        GC_ALLOCATE(BackboneDependentRotamerSet_O, me );
	me->_Phi = phi;
	me->_Psi = psi;
	return me;
    };

#else
    core::T_sp BackboneDependentRotamerSet_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
	this->_Phi = translate::from_object<int>::convert(env->lookup(_sym_phi));
	this->_Psi = translate::from_object<int>::convert(env->lookup(_sym_psi));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void BackboneDependentRotamerSet_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->attribute("phi",this->_Phi);
	node->attribute("psi",this->_Psi);
    }
#endif
    

    



    EXPOSE_CLASS(kinematics,BackboneDependentRotamerLibrary_O);
    void BackboneDependentRotamerLibrary_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<BackboneDependentRotamerLibrary_O>()
//	    .def_raw("__init__",&BackboneDependentRotamerLibrary_O::__init__,"(self &key phiStart phiStep phiCount psiStart psiStep psiCount)")
	    .def("backbone-dependent-rotamer-library-add-rotamer-set",&BackboneDependentRotamerLibrary_O::addRotamerSet)
	    .def("backbone-dependent-rotamer-library-rotamer-sets-as-cons",&BackboneDependentRotamerLibrary_O::rotamerSetsAsCons)
	    ;
//	    .def("backbone-dependent-rotamer-library-add-rotamer",&BackboneDependentRotamerLibrary_O::rotamerSetAddRotamer)
//	    .def("backbone-dependent-rotamer-library-as-cons",&BackboneDependentRotamerLibrary_O::rotamerSetAsCons)
//	    .def("backbone-dependent-rotamer-library-size",&BackboneDependentRotamerLibrary_O::size)
//	    .def("backbone-dependent-rotamer-library-get",&BackboneDependentRotamerLibrary_O::rotamerSetGet)
	;
    }

    void BackboneDependentRotamerLibrary_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,BackboneDependentRotamerLibrary,"","",_lisp)
//	    .def_raw("__init__",&"(self &key phiStart phiStep phiCount psiStart psiStep psiCount)")
	    .def("backbone-dependent-rotamer-library-add-rotamer-set",&BackboneDependentRotamerLibrary_O::addRotamerSet)
	    .def("backbone-dependent-rotamer-library-rotamer-sets-as-cons",&BackboneDependentRotamerLibrary_O::rotamerSetsAsCons)
	    ;
#endif
//	    .def("backbone-dependent-rotamer-library-add-rotamer",&BackboneDependentRotamerLibrary_O::rotamerSetAddRotamer)
//	    .def("backbone-dependent-rotamer-library-as-cons",&BackboneDependentRotamerLibrary_O::rotamerSetAsCons)
//	    .def("backbone-dependent-rotamer-library-size",&BackboneDependentRotamerLibrary_O::size)
//	    .def("backbone-dependent-rotamer-library-get",&BackboneDependentRotamerLibrary_O::rotamerSetGet)
	;
    }


#if INIT_TO_FACTORIES

#define ARGS_BackboneDependentRotamerLibrary_O_make "(phi-step phi-start phi-count psi-step psi-start psi-count)"
#define DECL_BackboneDependentRotamerLibrary_O_make ""
#define DOCS_BackboneDependentRotamerLibrary_O_make "make BackboneDependentRotamerLibrary"
  BackboneDependentRotamerLibrary_sp BackboneDependentRotamerLibrary_O::make(const int phiStep, const int phiStart, const int phiCount, const int psiStep, const int psiStart, const int psiCount)
    {_G();
	GC_ALLOCATE(BackboneDependentRotamerLibrary_O, me );
      me->_PhiStep = phiStep;
      me->_PhiStart = phiStart;
      me->_PhiCount = phiCount;
      me->_PsiStep = psiStep;
      me->_PsiStart = psiStart;
      me->_PsiCount = psiCount;
      ASSERTF(me->_PhiStart>=-180 && me->_PhiStart<180,BF("Illegal phiStart[%d] - must be in [-180,180)") % me->_PhiStart);
      ASSERTF(me->_PhiCount>0,BF("Illegal phiCount[%d]") % me->_PhiCount );
      ASSERTF(me->_PhiStep>0 && me->_PhiStep<=180,BF("Illegal phiStep[%d] - must be in [1,180)") % me->_PhiStep);
      ASSERTF(me->_PsiStart>=-180 && me->_PsiStart<180,BF("Illegal psiStart[%d] - must be in [-180,180)") % me->_PsiStart);
      ASSERTF(me->_PsiCount>0,BF("Illegal psiCount[%d]") % me->_PsiCount );
      ASSERTF(me->_PsiStep>0 && me->_PsiStep<=180,BF("Illegal psiStep[%d] - must be in [1,180)") % me->_PsiStep);
      me->_RotamerSets.resize(me->_PhiCount*me->_PsiCount);
      return me;
    };

#else
    core::T_sp BackboneDependentRotamerLibrary_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
	this->_PhiStep = translate::from_object<int>::convert(env->lookup(Pkg(),"phiStep"));
	this->_PhiStart = translate::from_object<int>::convert(env->lookup(Pkg(),"phiStart"));
	this->_PhiCount = translate::from_object<int>::convert(env->lookup(Pkg(),"phiCount"));
	this->_PsiStep = translate::from_object<int>::convert(env->lookup(Pkg(),"psiStep"));
	this->_PsiStart = translate::from_object<int>::convert(env->lookup(Pkg(),"psiStart"));
	this->_PsiCount = translate::from_object<int>::convert(env->lookup(Pkg(),"psiCount"));
	ASSERTF(this->_PhiStart>=-180 && this->_PhiStart<180,BF("Illegal phiStart[%d] - must be in [-180,180)") % this->_PhiStart);
	ASSERTF(this->_PhiCount>0,BF("Illegal phiCount[%d]") % this->_PhiCount );
	ASSERTF(this->_PhiStep>0 && this->_PhiStep<=180,BF("Illegal phiStep[%d] - must be in [1,180)") % this->_PhiStep);
	ASSERTF(this->_PsiStart>=-180 && this->_PsiStart<180,BF("Illegal psiStart[%d] - must be in [-180,180)") % this->_PsiStart);
	ASSERTF(this->_PsiCount>0,BF("Illegal psiCount[%d]") % this->_PsiCount );
	ASSERTF(this->_PsiStep>0 && this->_PsiStep<=180,BF("Illegal psiStep[%d] - must be in [1,180)") % this->_PsiStep);
	this->_RotamerSets.resize(this->_PhiCount*this->_PsiCount);
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void BackboneDependentRotamerLibrary_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->attribute("phiStart",this->_PhiStart);
	node->attribute("phiStep",this->_PhiStep);
	node->attribute("phiCount",this->_PhiCount);
	node->attribute("psiStart",this->_PsiStart);
	node->attribute("psiStep",this->_PsiStep);
	node->attribute("psiCount",this->_PsiCount);
	node->archiveVector0("rotamerSets",this->_RotamerSets);

    }
#endif


    void BackboneDependentRotamerLibrary_O::initialize()
    {_OF();
        this->Base::initialize();
    }


    bool BackboneDependentRotamerLibrary_O::validPhiPsi(int phi, int psi) const
    {_OF();
	int ophi = phi-this->_PhiStart;
	int mphi = ophi/this->_PhiStep;
	int rphi = ophi % this->_PhiStep;
	if ( rphi != 0 ) return false;
	if ( mphi > this->_PhiCount ) return false;
	int opsi = psi-this->_PsiStart;
	int mpsi = opsi/this->_PsiStep;
	int rpsi = opsi % this->_PsiStep;
	if ( rpsi != 0 ) return false;
	if ( mpsi > this->_PsiCount ) return false;
	return true;
    }

    int	BackboneDependentRotamerLibrary_O::_index(int phi, int psi, int& mphi, int& mpsi) const
    {_OF();
	int ophi = phi-this->_PhiStart;
	mphi = ophi/this->_PhiStep;
	int opsi = psi-this->_PsiStart;
	mpsi = opsi/this->_PsiStep;
	int idx = mphi*this->_PhiCount+mpsi;
	return idx;
    }
	

    void BackboneDependentRotamerLibrary_O::addRotamerSet(BackboneDependentRotamerSet_sp rotamerSet)
    {_OF();
	if ( !this->validPhiPsi(rotamerSet->_Phi,rotamerSet->_Psi) )
	{
	    int mphi;
	    int mpsi;
	    this->_index(rotamerSet->_Phi,rotamerSet->_Psi,mphi,mpsi);
	    SIMPLE_ERROR(BF("Illegal phi[%d] or psi[%d] mphi[%d] mpsi[%d] phiCount[%d] psiCount[%d] for BackboneDependentRotamerLibrary")
			       % rotamerSet->_Phi % rotamerSet->_Psi
			       % mphi % mpsi
			       % this->_PhiCount % this->_PsiCount
		      );
	}
	int mphi, mpsi;
	int idx = this->_index(rotamerSet->_Phi,rotamerSet->_Psi,mphi,mpsi);
	this->_RotamerSets[idx] = rotamerSet;
    }

    core::Cons_sp BackboneDependentRotamerLibrary_O::rotamerSetsAsCons() const
    {_OF();
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>());
	core::Cons_sp cur = first;
	for ( int i=0; i<(int)this->_RotamerSets.size(); i++ )
	{
	    core::Cons_sp one = core::Cons_O::create(this->_RotamerSets[i]);
	    cur->setCdr(one);
	    cur = one;
	}
	return first->cdr();
    }

}; /* kinematics */
