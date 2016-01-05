       
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/plug.h>
#include <cando/chem/candoDatabase.h>
#include <clasp/core/cons.h>
#include <cando/chem/coupling.h>
#include <clasp/core/environment.h>
#include <clasp/core/numbers.h>
#include <cando/chem/constitution.h>
#include <clasp/core/wrappers.h>


namespace chem
{



#if 0
    RingClosingMate_sp RingClosingMate_O::create(CandoDatabase_sp db)
    {
	GC_ALLOCATE(RingClosingMate_O, mate );
//    mate->setCandoDatabase(db);
	return mate;
    }
#endif


    core::Symbol_sp RingClosingMate_O::getName() const
    {
	return _lisp->intern("UnnamedRingClosingMate");
    };


//
// Constructor
//


//
// Copy Constructor
//
    RingClosingMate_O::RingClosingMate_O(const RingClosingMate_O& ss) : Base(ss)
    {
    }





    string RingClosingMate_O::descriptionOfContents() const
    {
	stringstream ss, shex;
	return ss.str();
    }

#ifdef XML_ARCHIVE
    void	RingClosingMate_O::archiveBase(core::ArchiveP node)
    {_G();
	this->Base::archiveBase(node);
    }
#endif

/*
  __BEGIN_DOC(classes.RingClosingMate)
  \inheritsFrom{EntityNameSetWithCap}
  A EntityNameSet that keeps track of a capping monomer that is used to cap training oligomers when they are being defined.
  The capping monomer is supposed to be small and best represent the other members of the RingClosingMate.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.RingClosingMate.!class.RingClosingMate)
  \requiredKeyed{groupNames:}{List::groupNames}

  Initialize a RingClosingMate object.
  __END_DOC
*/

    void Mate_O::initialize()
    {
	// nothing 
    }

    core::Symbol_sp Mate_O::getName() const
    {_G();
	return _Nil<core::Symbol_O>();
    }


#if 0
    Mate_sp Mate_O::create(CandoDatabase_sp db)
    {
	GC_ALLOCATE(Mate_O, mate );
	mate->setCandoDatabase(db);
	return mate;
    }
#endif

//
// Constructor
//


//
// Copy Constructor
//
    Mate_O::Mate_O(const Mate_O& ss) :	Base(ss)
    {
	this->_Cap = ss._Cap;
    }





    string Mate_O::descriptionOfContents() const
    {
	stringstream ss, shex;
	ss << "( Cap "<< this->_Cap<<") ";
	return ss.str();
    }

#ifdef XML_ARCHIVE
    void	Mate_O::archiveBase(core::ArchiveP node)
    {_G();
	this->Base::archiveBase(node);
	node->attribute( "cap", this->_Cap );
    }
#endif

/*
  __BEGIN_DOC(classes.Mate)
  \inheritsFrom{EntityNameSetWithCap}
  A EntityNameSet that keeps track of a capping monomer that is used to cap training oligomers when they are being defined.
  The capping monomer is supposed to be small and best represent the other members of the Mate.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.Mate.!class.Mate)
  \requiredKeyed{cap:}{Text::capName}
  \requiredKeyed{groupNames:}{List::groupNames}

  Initialize a Mate object.
  __END_DOC
*/

#if INIT_TO_FACTORIES

#define ARGS_Mate_O_make "(cap)"
#define DECL_Mate_O_make ""
#define DOCS_Mate_O_make "make Mate"
  Mate_sp Mate_O::make(core::Symbol_sp cap)
  {_G();
      GC_ALLOCATE(Mate_O, me );
    me->_Cap = cap;
    return me;
  };

#else

    core::T_sp Mate_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	this->Base::__init__(exec,args,env,lisp);
	this->_Cap = translate::from_object<core::Symbol_O>::convert(env->lookup(ChemPkg,"cap"));
	return _Nil<core::T_O>();
    }

#endif





    void	Plug_O::initialize()
    {
	this->Base::initialize();
	this->_WeakConstitution = _Nil<Constitution_O>();
	this->_Name = _Nil<core::Symbol_O>();
    }


/*
  __BEGIN_DOC(classes.Plug)
  Defines one or two atoms of this monomer that can be plugged into, a plug name and a
  collection of Mate objects that can act as mates for this plug.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.Plug.!class.Plug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a Plug object. Plugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  Outgoing plugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  __END_DOC
*/

#if INIT_TO_FACTORIES

#define ARGS_Plug_O_make "(name)"
#define DECL_Plug_O_make ""
#define DOCS_Plug_O_make "make Plug"
  Plug_sp Plug_O::make(core::Symbol_sp name)
  {_G();
      GC_ALLOCATE(Plug_O, me );
    me->_Name = name;
    return me;
  };

#else

    core::T_sp Plug_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	this->_Name = translate::from_object<core::Symbol_O>::convert(env->lookup(ChemPkg,"name"));
	return _Nil<core::T_O>();
    }

#endif

    Plug_O::Plug_O(const Plug_O& p) : core::CxxObject_O(p)
    {_G();
	Mate_sp	rn;
	LOG(BF("Original %s") % p.description().c_str()  );
	this->_WeakConstitution = p._WeakConstitution;
	this->_Name = p._Name;
	LOG(BF("Copy constructed %s") % this->description().c_str()  );
    }



    string	Plug_O::descriptionOfContents() const
    {
	stringstream ss;
	ss << this->Base::descriptionOfContents() << " ";
	ss << "( Name " << this->_Name->__repr__() << " )";
	return ss.str();
    }




    void	Plug_O::setConstitution(Constitution_sp con)
    {
	this->_WeakConstitution = con->sharedThis<Constitution_O>();
    }



CL_LISPIFY_NAME("getConstitution");
CL_DEFMETHOD     Constitution_sp	Plug_O::getConstitution()
    {_OF();
	ASSERTNOTNULL(this->_WeakConstitution);
	if ( this->_WeakConstitution.nilp() )
	{
	    // Look for a constitution in the owners
	    IMPLEMENT_MEF(BF("Handle ownerWithClass"));
//	    this->_WeakConstitution = this->ownerWithClass<Constitution_O>();
	}
	return this->_WeakConstitution;
    }





CL_LISPIFY_NAME("otherSidePlugName");
CL_DEFMETHOD     core::Symbol_sp Plug_O::otherSidePlugName()
    {
	return DirectionalCoupling_O::otherPlugName(this->_Name);
    }


#ifdef XML_ARCHIVE
    void	Plug_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	if ( node->saving() ) this->getConstitution();
	node->archiveWeakPointer("constitution",this->_WeakConstitution);
#if PRODUCTION_CODE   // FIXME use "name" only and remove the test for "_key"
	node->attribute("name",this->_Name);
#else
	if ( node->loading() )
	{
	    if ( node->hasAttribute("name") )
	    {
		node->attribute("name",this->_Name);
	    } else
	    {
		node->attribute("_key",this->_Name);
	    }
	} else
	{
	    node->attribute("name",this->_Name);
	}
#endif
    }
#endif


    gctools::Vec0<Mate_sp>	_EmptyMateList;

    gctools::Vec0<Mate_sp>::iterator	Plug_O::begin_Mates()
    {
	return _EmptyMateList.end();
    }

    gctools::Vec0<Mate_sp>::iterator	Plug_O::end_Mates()
    {
	return _EmptyMateList.end();
    }














/*
  __BEGIN_DOC(classes.PlugWithMates)
  Defines one or two atoms of this monomer that can be plugged into, a plug name and a
  collection of Mate objects that can act as mates for this plug.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.PlugWithMates.!class.PlugWithMates)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a PlugWithMates object. PlugWithMatess can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  __END_DOC
*/
    void PlugWithMates_O::initialize()
    {
	this->Base::initialize();
    }

#if INIT_TO_FACTORIES

#define ARGS_PlugWithMates_O_make "(bond0 bond1 mates)"
#define DECL_PlugWithMates_O_make ""
#define DOCS_PlugWithMates_O_make "make PlugWithMates"
  PlugWithMates_sp PlugWithMates_O::make(core::Symbol_sp bond0, core::Symbol_sp bond1, core::List_sp mates)
  {_G();
      GC_ALLOCATE(PlugWithMates_O, me );
      me->_B0 = bond0;
      me->_B1 = bond1;
      core::fillVec0(mates,me->_Mates);
      return me;
  };

#else

    core::T_sp PlugWithMates_O::__init__(core::Function_sp exec, core::List_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	this->Base::__init__(exec,args,env,lisp);
	this->_B0 = translate::from_object<string>::convert(env->lookup(ChemPkg,"bond0"));
	this->_B1 = translate::from_object<string>::convert(env->lookup(ChemPkg,"bond1"));
	core::List_sp mates = translate::from_object<core::List_O>::convert(env->lookup(ChemPkg,"mates"));
        mates->fillVec0(this->_Mates);
	return _Nil<core::T_O>();
    }
#endif


    PlugWithMates_O::PlugWithMates_O(const PlugWithMates_O& p) : Plug_O(p)
    {_G();
	LOG(BF("Original %s") % p.description().c_str()  );
	this->_B0 = p._B0;
	this->_B1 = p._B1;
	this->_Mates.clear();
        gctools::Vec0<Mate_sp>::const_iterator vi;
	LOG(BF("Copying %d mates") % p._Mates.size()  );
	for ( vi=p._Mates.begin(); vi!=p._Mates.end(); vi++ ) {
	    LOG(BF("Copied mate") );
	    GC_COPY(Mate_O, rn, *vi->get() ); // = RP_Copy<Mate_O>(*vi);
	    this->_Mates.push_back(rn);
	}
	LOG(BF("Copy constructed %s") % this->description().c_str()  );
    }

    string	PlugWithMates_O::descriptionOfContents() const
    {
	stringstream ss;
	ss << this->Base::descriptionOfContents() << " ";
	ss << " ( Mates ";
	ss << "#" << this->_Mates.size() << " ";
        gctools::Vec0<Mate_sp>::const_iterator vi;
	for ( vi=this->_Mates.begin(); vi!=this->_Mates.end(); vi++ )
	{
	    ss << (*vi)->description() << ", ";
	}
	ss << " )";
	return ss.str();
    }


    string PlugWithMates_O::__repr__() const
    {
	stringstream ss;
	ss << "( " << this->className() << " '"<< this->getName()->__repr__() << " :mates ";
	ss << "(list ";
	for ( Mates::const_iterator it=this->_Mates.begin(); it!=this->_Mates.end(); it++ )
	{
	    ss << (*it)->__repr__() << " ";
	}
	return ss.str();
    }

#ifdef XML_ARCHIVE
    void	PlugWithMates_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attributeIfNotDefault<string>("b0",this->_B0,"");
	node->attributeIfNotDefault<string>("b1",this->_B1,"");
	node->archiveVector0("mates",this->_Mates);
    }
#endif



    bool	PlugWithMates_O::recognizesMateNameOrPdb(core::Symbol_sp name)
    {_G();
        gctools::Vec0<Mate_sp>::iterator	mi;
	for ( mi=this->_Mates.begin(); mi!=this->_Mates.end(); mi++ ) {
	    if ( (*mi)->recognizesNameOrPdb(name) ) return true;
	}
	return false;
    }

    void OutPlug_O::initialize()
    {
	this->Base::initialize();
    }



/*
  __BEGIN_DOC(classes.OutPlug.!class.OutPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a OutPlug object. OutPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  Outgoing plugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  __END_DOC
*/

#if INIT_TO_FACTORIES

#define ARGS_OutPlug_O_make "(stub_pivot_atom)"
#define DECL_OutPlug_O_make ""
#define DOCS_OutPlug_O_make "make OutPlug"
  OutPlug_sp OutPlug_O::make(MatterName stubPivotAtom)
  {_G();
      GC_ALLOCATE(OutPlug_O, me );
    me->_StubPivotAtom = stubPivotAtom;
    return me;
  };

#else

    core::T_sp OutPlug_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	this->Base::__init__(exec,args,env,lisp);
	this->_StubPivotAtom = translate::from_object<string>::convert(env->lookup(ChemPkg,"stubPivotAtom"));
	return _Nil<core::T_O>();
    }

#endif

    OutPlug_O::OutPlug_O(const OutPlug_O& p) : PlugWithMates_O(p)
    {_G();
	Mate_sp	rn;
	LOG(BF("Original %s") % p.description().c_str()  );
	IMPLEMENT_ME(); // Should the next copy be shallow or deep?
	LOG(BF("Copy constructed %s") % this->description().c_str()  );
    }



    string	OutPlug_O::descriptionOfContents() const
    {
	stringstream ss, shex;
	ss << this->Base::descriptionOfContents() << " ";
	return ss.str();
    }



#ifdef XML_ARCHIVE
    void	OutPlug_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute("stubPivotAtom",this->_StubPivotAtom);
    }
#endif



    void InPlug_O::initialize()
    {
	this->Base::initialize();
    }


/*
  __BEGIN_DOC(classes.InPlug.!class.InPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a InPlug object. InPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  __END_DOC
*/
    InPlug_O::InPlug_O(const InPlug_O& p) : PlugWithMates_O(p)
    {_G();
	Mate_sp	rn;
	LOG(BF("Original %s") % p.description().c_str()  );
    }



    string	InPlug_O::descriptionOfContents() const
    {
	return this->Base::descriptionOfContents();
    }







    void JumpPlug_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<JumpPlug_O>()
//	    .def_raw("core:__init__",&JumpPlug_O::__init__,"(self &key name jumpAtomName)") // inherit from Plug_O.initArgs("(self &key name)")
	    .def("rootAtomName",&JumpPlug_O::rootAtomName)
	    ;
    }

    void JumpPlug_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,JumpPlug,"","",_lisp)
	    .def("rootAtomName",&JumpPlug_O::rootAtomName)
	    ;
#endif
    }





    void JumpPlug_O::initialize()
    {
	this->Base::initialize();
	this->_JumpAtomName = _Nil<core::Symbol_O>();
    }

/*
  __BEGIN_DOC(classes.JumpPlug.!class.JumpPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{originFrame:}{Frame::originFrame}

  Initialize a JumpPlug object. JumpPlugs don't make bonds don't have mates but they do have an origin frame that is within the topology.
  __END_DOC
*/

#if INIT_TO_FACTORIES

#define ARGS_JumpPlug_O_make "(jump_atom_name)"
#define DECL_JumpPlug_O_make ""
#define DOCS_JumpPlug_O_make "make JumpPlug"
  JumpPlug_sp JumpPlug_O::make(MatterName jumpAtomName)
  {_G();
      GC_ALLOCATE(JumpPlug_O, me );
    me->_JumpAtomName = jumpAtomName;
    return me;
  };

#else

    core::T_sp JumpPlug_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	this->Base::__init__(exec,args,env,lisp);
	this->_JumpAtomName = translate::from_object<string>::convert(env->lookup(ChemPkg,"jumpAtomName"));
	return _Nil<core::T_O>();
    }

#endif


    JumpPlug_O::JumpPlug_O(const JumpPlug_O& p) : Plug_O(p)
    {_G();
	Mate_sp	rn;
	this->_JumpAtomName = p._JumpAtomName;
	LOG(BF("Jumpal %s") % p.description().c_str()  );
    }



    string	JumpPlug_O::descriptionOfContents() const
    {
	stringstream ss;
	ss << this->Base::descriptionOfContents() << " ";
	ss << "JumpAtomName[" << this->_JumpAtomName << "]";
	return ss.str();
    }


#ifdef XML_ARCHIVE
    void	JumpPlug_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute("jumpAtom",this->_JumpAtomName);
    }
#endif














    void RingClosingPlug_O::initialize()
    {
	this->Base::initialize();
    }


/*
  __BEGIN_DOC(classes.RingClosingPlug.!class.RingClosingPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}
  \optionalKeyed{exportFrame:}{Frame::exportFrame}
  \requiredKeyed{ringClosingMates:}{Cons::listOfRingClosingMates}

  Initialize a RingClosingPlug object. RingClosingPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  RingClosingPlugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  There is an additional list of ringClosingMates that can be attached to this plug without building the mates.
  __END_DOC
*/

#if INIT_TO_FACTORIES

#define ARGS_RingClosingPlug_O_make "(ring_closing_mates)"
#define DECL_RingClosingPlug_O_make ""
#define DOCS_RingClosingPlug_O_make "make RingClosingPlug"
  RingClosingPlug_sp RingClosingPlug_O::make(core::List_sp ringClosingMates)
  {_G();
      GC_ALLOCATE(RingClosingPlug_O, me );
      core::fillVec0(gc::As<core::Cons_sp>(ringClosingMates),me->_RingClosingMates);
    return me;
  };

#else

    core::T_sp RingClosingPlug_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	this->Base::__init__(exec,args,env,lisp);
	core::Cons_sp mates = translate::from_object<core::Cons_sp>::convert(env->lookup(ChemPkg,"ringClosingMates"));
	this->_RingClosingMates.fillFromCons(mates);
	return _Nil<core::T_O>();
    }

#endif

    RingClosingPlug_O::RingClosingPlug_O(const RingClosingPlug_O& p) : RingClosingPlug_O::Base(p)
    {_G();
	LOG(BF("Original %s") % p.description().c_str()  );
        gctools::Vec0<RingClosingMate_sp>::const_iterator vi;
	LOG(BF("Copying %d ring closing mates") % p._Mates.size()  );
	for ( vi=p._RingClosingMates.begin(); vi!=p._RingClosingMates.end(); vi++ ) {
	    LOG(BF("Copied mate") );
	    GC_COPY(RingClosingMate_O, rn, *vi->get() ); // = RP_Copy<RingClosingMate_O>(*vi);
	    this->_RingClosingMates.push_back(rn);
	}
	LOG(BF("Copy constructed %s") % this->description().c_str()  );
    }



    string	RingClosingPlug_O::descriptionOfContents() const
    {
	stringstream ss, shex;
	ss << this->Base::descriptionOfContents() << " ";
	ss << " ( RingClosingMates ";
	ss << "#" << this->_RingClosingMates.size() << " ";
        gctools::Vec0<RingClosingMate_sp>::const_iterator vi;
	for ( vi=this->_RingClosingMates.begin(); vi!=this->_RingClosingMates.end(); vi++ ) {
	    ss << (*vi)->description() << ", ";
	}
	ss << " )";
	return ss.str();
    }


#ifdef XML_ARCHIVE
    void	RingClosingPlug_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->archiveVector0("ringCLosingMates",this->_RingClosingMates);
    }
#endif








    void RingClosingMate_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<RingClosingMate_O>()
//	    .def_raw("core:__init__",&RingClosingMate_O::__init__,"(self &key cap)")
	    ;
    }

    void RingClosingMate_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,RingClosingMate,"","",_lisp)
	    ;
#endif
    }





    void Mate_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Mate_O>()
//	    .def_raw("core:__init__",&Mate_O::__init__,"(self &key entityNames optional cap)")
	    .def("setCap",&Mate_O::setCap)
	    .def("getCap",&Mate_O::getCap)
	    ;
    }

    void Mate_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Mate,"","",_lisp)
	    .def("setCap",&Mate_O::setCap)
	    .def("getCap",&Mate_O::getCap)
	    ;
#endif
    }



    string Mate_O::__repr__() const
    {
	stringstream ss;
	ss << "( " << this->className();
	{
	    if ( this->getName().notnilp() )
	    {
		ss << " :name " << this->getName()->__repr__();
	    }
	}
	ss << _rep_(this->_EntityNames);
	ss << " )";
	return ss.str();

    }


    void Plug_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Plug_O>()
//	    .def_raw("core:__init__",&Plug_O::__init__,"(self &key name)")
	    .def("getConstitution", &Plug_O::getConstitution)
	    .def("getName", &Plug_O::getName,"","","",false)
	    .def("setName", &Plug_O::setName,"","","",false)
	    .def("otherSidePlugName", &Plug_O::otherSidePlugName)
	    .def("getPlug", &Plug_O::getPlug)
	    .def("getIsIn", &Plug_O::getIsIn)
	    .def("getIsRingClosing", &Plug_O::getIsRingClosing)
	    .def("matesAsCons",&Plug_O::matesAsCons)
	    ;
    }

    void Plug_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Plug,"","",_lisp)
	    .def("getConstitution", &Plug_O::getConstitution)
	    .def("getName", &Plug_O::getName)
	    .def("setName", &Plug_O::setName)
	    .def("otherSidePlugName", &Plug_O::otherSidePlugName)
	    .def("getPlug", &Plug_O::getPlug)
	    .def("getIsIn", &Plug_O::getIsIn)
	    .def("getIsRingClosing", &Plug_O::getIsRingClosing)
	    .def("matesAsCons",&Plug_O::matesAsCons)
	    ;
#endif
    }


    void PlugWithMates_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<PlugWithMates_O>()
//	    .def_raw("core:__init__",&PlugWithMates_O::__init__,"(self &key name bond0 (bond1 \"\") mates)")
	    .def("getB0",&PlugWithMates_O::getB0)
	    .def("getB1",&PlugWithMates_O::getB1)
	    .def("hasB1",&PlugWithMates_O::hasB1)
	    ;
    }

    void PlugWithMates_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,PlugWithMates,"","",_lisp)
	    .def("getB0",&PlugWithMates_O::getB0)
	    .def("getB1",&PlugWithMates_O::getB1)
	    .def("hasB1",&PlugWithMates_O::hasB1)
	    ;
#endif
    }





    void OutPlug_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<OutPlug_O>()
//	    .def_raw("core:__init__",&OutPlug_O::__init__,"(self &key name bond0 (bond1 \"\") mates (stubPivotAtom \"\"))") // inherits from PlugWithMates_O.initArgs("(self &key name bond0 (bond1 \"\") mates)")
	    .def("hasStubPivotAtom",&OutPlug_O::hasStubPivotAtom)
	    .def("getStubPivotAtom",&OutPlug_O::getStubPivotAtom)
	    ;
    }

    void OutPlug_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,OutPlug,"","",_lisp)
	    .def("hasStubPivotAtom",&OutPlug_O::hasStubPivotAtom)
	    .def("getStubPivotAtom",&OutPlug_O::getStubPivotAtom)
	    ;
#endif
    }




    void InPlug_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<InPlug_O>()
//	    .def_raw("core:__init__",&InPlug_O::__init__,"(self &key name bond0 (bond1 \"\") mates)") // inherits from PlugWithMates_O.initArgs("(self &key name bond0 (bond1 \"\") mates)")
	    .def("rootAtomName",&InPlug_O::rootAtomName)
	    ;
    }

    void InPlug_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,InPlug,"","",_lisp)
	    .def("rootAtomName",&InPlug_O::rootAtomName)
	    ;
#endif
    }



    void RingClosingPlug_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<RingClosingPlug_O>()
//	    .def_raw("core:__init__",&RingClosingPlug_O::__init__,"(self &key name bond0 (bond1 \"\") mates (stubPivotAtom \"\"))")
	    .def("ringClosingMatesAsCons", &RingClosingPlug_O::ringClosingMatesAsCons)
	    ;
    }

    void RingClosingPlug_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,RingClosingPlug,"","",_lisp)
	    .def("ringClosingMatesAsCons", &RingClosingPlug_O::ringClosingMatesAsCons)
	    ;
#endif
    }




CL_LISPIFY_NAME("ringClosingMatesAsCons");
CL_DEFMETHOD     core::List_sp RingClosingPlug_O::ringClosingMatesAsCons()
    {_OF();
	return core::Cons_O::createFromVec0(this->_RingClosingMates);
    }


    bool RingClosingPlug_O::recognizesRingClosingMate(core::Symbol_sp mateName)
    {_OF();
        gctools::Vec0<RingClosingMate_sp>::iterator it;
	for ( it=this->_RingClosingMates.begin(); it!=this->_RingClosingMates.end(); it++ )
	{
	    if ( (*it)->recognizesMonomerName(mateName)) return true;
	}
	return false;
    }

    EXPOSE_CLASS(chem, Mate_O);
    EXPOSE_CLASS(chem, RingClosingMate_O);
    EXPOSE_CLASS(chem, Plug_O);
    EXPOSE_CLASS(chem, PlugWithMates_O);
    EXPOSE_CLASS(chem, OutPlug_O);
    EXPOSE_CLASS(chem, InPlug_O);
    EXPOSE_CLASS(chem, JumpPlug_O);
    EXPOSE_CLASS(chem, RingClosingPlug_O);


}; //namespace chem


