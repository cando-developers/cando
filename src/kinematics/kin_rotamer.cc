/*
    File: kin_rotamer.cc
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

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/rotamer.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//





#if INIT_TO_FACTORIES

CL_LISPIFY_NAME(make-Rotamer);
Rotamer_sp Rotamer_O::make(core::List_sp dihedrals, core::List_sp sigmas, core::List_sp indices, const double probability, const int count)
    {
      auto  me  = gctools::GC<Rotamer_O>::allocate_with_default_constructor();
      me->_Count = count;
      me->_Probability = probability;
      for ( ; dihedrals.notnilp(); dihedrals = oCdr(dihedrals),
	      sigmas = oCdr(sigmas),
	      indices = oCdr(indices) )
	{
	  me->_Dihedrals.push_back(oCar(dihedrals).as<core::DoubleFloat_O>()->get());
	  me->_Sigmas.push_back(oCar(sigmas).as<core::DoubleFloat_O>()->get());
	  if ( oCar(indices).notnilp() )
	    {
	      me->_Indices.push_back(oCar(indices).as<core::Fixnum_I>().unsafe_fixnum());
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
    {
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
	    this->_Dihedrals.push_back(dihedrals->car<core::DoubleFloat_O>());
	    this->_Sigmas.push_back(sigmas->car<core::DoubleFloat_O>());
	    if ( indices->ocar().notnilp() )
	    {
		this->_Indices.push_back(indices->car<core::Fixnum_O>());
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








#if INIT_TO_FACTORIES

CL_LISPIFY_NAME(make-ChiDihedral);
  RotamerSetBase_sp RotamerSetBase_O::make()
    {
      IMPLEMENT_ME();
    };

#else
    core::T_sp RotamerSetBase_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
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
	this->_Rotamers.push_back(rotamer);
	return rotamer;
    }

    Rotamer_sp RotamerSetBase_O::get(int i) const
    {_OF();
	ASSERTF(i>=0 && i<(int)this->_Rotamers.size(),BF("Illegal index for rotamer[%d]") % i);
	return this->_Rotamers[i];
    }
    

core::List_sp RotamerSetBase_O::asList() const
{
  core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
  core::Cons_sp cur = first;
  for ( int i=0; i<(int)this->_Rotamers.size(); i++ )
  {
    core::Cons_sp one = core::Cons_O::create(this->_Rotamers[i],nil<core::T_O>());
    cur->setCdr(one);
    cur = one;
  }
  return oCdr(first);
}




















#if INIT_TO_FACTORIES

CL_LISPIFY_NAME(make-RotamerSet)
RotamerSet_sp RotamerSet_O::make()
    {
      IMPLEMENT_ME();
    };

#else
    core::T_sp RotamerSet_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
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











    void BackboneDependentRotamerSet_O::initialize()
    {
	this->Base::initialize();
    }

#if INIT_TO_FACTORIES

CL_LISPIFY_NAME(make-BackboneDependentRotamerSet);
BackboneDependentRotamerSet_sp BackboneDependentRotamerSet_O::make(const int phi, const int psi)
    {
      auto  me  = gctools::GC<BackboneDependentRotamerSet_O>::allocate_with_default_constructor();
	me->_Phi = phi;
	me->_Psi = psi;
	return me;
    };

#else
    core::T_sp BackboneDependentRotamerSet_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
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
    

    






#if INIT_TO_FACTORIES

CL_LISPIFY_NAME(make-BackboneDependentRotamerLibrary);
BackboneDependentRotamerLibrary_sp BackboneDependentRotamerLibrary_O::make(const int phiStep, const int phiStart, const int phiCount, const int psiStep, const int psiStart, const int psiCount)
    {
      auto  me  = gctools::GC<BackboneDependentRotamerLibrary_O>::allocate_with_default_constructor();
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
    {
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

    core::List_sp BackboneDependentRotamerLibrary_O::rotamerSetsAsList() const
    {
      core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
      core::Cons_sp cur = first;
      for ( int i=0; i<(int)this->_RotamerSets.size(); i++ )
      {
        core::Cons_sp one = core::Cons_O::create(this->_RotamerSets[i],nil<core::T_O>());
        cur->setCdr(one);
        cur = one;
      }
      return oCdr(first);
    }

}; /* kinematics */
