#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/conformation.h>
#include <cando/kinematics/foldTree.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/atomTree.h>
#include <cando/chem/oligomer.h>
#include <clasp/core/wrappers.h>


namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS_AND_GLOBALS(kinematics,Conformation_O);

    void Conformation_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();
    }


    void Conformation_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Conformation_O>()
//	    .def_raw("__init__",&Conformation_O::__init__,"(self &key oligomers)")
	    .def("resizeMolecules",&Conformation_O::resizeMolecules)
	    .def("buildMoleculeUsingOligomer",&Conformation_O::buildMoleculeUsingOligomer)
	    .def("getFoldTree",&Conformation_O::getFoldTree)
	    .def("getAtomTree",&Conformation_O::getAtomTree)
	;
    }

    void Conformation_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,Conformation,"","",_lisp)
//	    .def_raw("__init__",&"(self &key oligomers)")
	    .def("resizeMolecules",&Conformation_O::resizeMolecules)
	    .def("buildMoleculeUsingOligomer",&Conformation_O::buildMoleculeUsingOligomer)
	    .def("getFoldTree",&Conformation_O::getFoldTree)
	    .def("getAtomTree",&Conformation_O::getAtomTree)
	;
#endif
    }


#if INIT_TO_FACTORIES

#define ARGS_Conformation_O_make "(id oligomers)"
#define DECL_Conformation_O_make ""
#define DOCS_Conformation_O_make "make Conformation"
  Conformation_sp Conformation_O::make(const core::Cons_sp& oligomers)
    {_G();
	GC_ALLOCATE(Conformation_O, me );
      me->resizeMolecules(oligomers->length());
      int moleculeId = 0;
      for ( core::Cons_sp cur=oligomers; cur.notnilp(); cur = cur->cdr(),moleculeId++ )
      {
	  chem::Oligomer_sp oligomer = cur->car<chem::Oligomer_O>();
	  me->buildMoleculeUsingOligomer(moleculeId,oligomer);
      }
      return me;
    };

#else
    core::T_sp Conformation_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
      this->Base::__init__(exec,args,env,lisp);
      core::Cons_sp oligomerCons = translate::from_object<core::Cons_O>::convert(env->lookup(Pkg(),"oligomers"));
      this->resizeMolecules(oligomerCons->length());
      int moleculeId = 0;
      for ( core::Cons_sp cur=oligomerCons; cur.notnilp(); cur = cur->cdr(),moleculeId++ )
      {
	  chem::Oligomer_sp oligomer = cur->car<chem::Oligomer_O>();
	  this->buildMoleculeUsingOligomer(moleculeId,oligomer);
      }
      return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void Conformation_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif

    Conformation_sp Conformation_O::create(chem::Oligomer_sp oligomer)
    {
	GC_ALLOCATE(Conformation_O, conf );
	conf->resizeMolecules(1);
	conf->buildMoleculeUsingOligomer(0,oligomer);
	return conf;
    }


    void Conformation_O::initialize()
    {_OF();
        this->Base::initialize();
	this->_FoldTree = FoldTree_O::create();
	this->_AtomTree = AtomTree_O::create();
    }


    int Conformation_O::numberOfMolecules() const
    {_OF();
	return this->_FoldTree->numberOfChains();
    }

    void Conformation_O::resizeMolecules(int numMolecules)
    {_OF();
	ASSERTF(numMolecules > 0, BF("You must allocate more than 0 molecules"));
	int oldNumberOfMolecules = this->_FoldTree->numberOfChains();
	this->_FoldTree->resizeChains(numMolecules);
	this->_AtomTree->resizeMolecules(numMolecules);
	this->notify(_sym_resizeMoleculesEvent,
		     core::Cons_O::createList(kw::_sym_newValue,
					      core::clasp_make_fixnum(numMolecules),
					      kw::_sym_oldValue,
					      core::clasp_make_fixnum(oldNumberOfMolecules)));
    }


    void Conformation_O::buildMoleculeUsingOligomer(int moleculeId, chem::Oligomer_sp oligomer)
    {_OF();
	ASSERTF(moleculeId<this->numberOfMolecules(),
		BF("Illegal moleculeId[%d]") % moleculeId );
	ChainNode_sp chainNode =
	    this->_FoldTree->buildChainUsingOligomer(moleculeId,oligomer);
	LOG(BF("Built FoldTree--->\n%s")%chainNode->_RootMonomerNode->asString());
	this->_AtomTree->buildMoleculeUsingChainNode(moleculeId,chainNode,oligomer);
	this->notify(_sym_buildMoleculeUsingOligomerEvent,
		     core::Cons_O::create(core::clasp_make_fixnum(moleculeId),_lisp));
    }


}; /* kinematics */
