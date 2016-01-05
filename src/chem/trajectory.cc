#define	DEBUG_LEVEL_FULL

#include <clasp/core/hashTableEq.h>
#include <cando/chem/trajectory.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/loop.h>
#include <cando/geom/coordinateArray.h>

// last include is wrappers.h
#include <clasp/core/wrappers.h>


namespace chem {


    EXPOSE_CLASS(chem,TrajectoryFrame_O);
    EXPOSE_CLASS(chem,Trajectory_O);




void TrajectoryFrame_O::exposeCando(core::Lisp_sp e)
{
    core::class_<TrajectoryFrame_O>()
    ;
}

    void TrajectoryFrame_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,TrajectoryFrame,"","",_lisp)
    ;
#endif //]
}

#if INIT_TO_FACTORIES

#define ARGS_TrajectoryFrame_O_make "()"
#define DECL_TrajectoryFrame_O_make ""
#define DOCS_TrajectoryFrame_O_make "make TrajectoryFrame"
  TrajectoryFrame_sp TrajectoryFrame_O::make()
  {_G();
    IMPLEMENT_ME();
  };

#else

    core::T_sp 	TrajectoryFrame_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{
    IMPLEMENT_ME();
    	// your stuff here
}

#endif

void	TrajectoryFrame_O::initialize()
{
    this->Base::initialize();
}

#ifdef XML_ARCHIVE
    void	TrajectoryFrame_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("coords",this->_Coordinates);
}
#endif


    void TrajectoryFrame_O::fillFromMatter(gctools::Vec0<Atom_sp>& atomList )
{
    this->_Coordinates = geom::CoordinateArray_O::create(atomList.size());
    geom::CoordinateArray_O::iterator ci;
    gctools::Vec0<Atom_sp>::const_iterator ai;
    for ( ai=atomList.begin(), ci = this->_Coordinates->begin(); ai!=atomList.end(); ai++, ci++ )
    {
	(*ci) = (*ai)->getPosition();
    }
}



    void	TrajectoryFrame_O::applyToMatter(gctools::Vec0<Atom_sp>& atomList )
{_OF();
    ASSERT_eq(this->_Coordinates->size(),atomList.size());
    geom::CoordinateArray_O::iterator ci;
    gctools::Vec0<Atom_sp>::iterator ai;
    for ( ai=atomList.begin(), ci = this->_Coordinates->begin(); ai!=atomList.end(); ai++, ci++ )
    {
	(*ai)->setPosition(*ci);
    }
}


void Trajectory_O::exposeCando(core::Lisp_sp e)
{
    Trajectory_O::expose(core::class_<Trajectory_O>());
}

    void Trajectory_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,Trajectory,"","",_lisp)
	    ;
#endif //]
}

#if INIT_TO_FACTORIES

#define ARGS_Trajectory_O_make "()"
#define DECL_Trajectory_O_make ""
#define DOCS_Trajectory_O_make "make Trajectory"
  Trajectory_sp Trajectory_O::make()
  {_G();
    IMPLEMENT_ME();
  };

#else

    core::T_sp 	Trajectory_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{
    IMPLEMENT_ME();
#if 0
    this->Base::oldLispInitialize(kargs,env);
    Matter_sp matter = kargs->getAndRemove("matter").as<Matter_O>();
    this->_setupAtomList(matter);
#endif
}

#endif

void	Trajectory_O::initialize()
{
    this->Base::initialize();
    this->_Matter = _Nil<Matter_O>();
    this->_Namespace = core::HashTableEq_O::create_default();
}


void	Trajectory_O::_setupAtomList(Matter_sp matter)
{
    this->_Matter = matter;
    this->_AtomList.clear();
    Loop lAtoms;
    lAtoms.loopTopGoal(this->_Matter,ATOMS);
    while ( lAtoms.advance() )
    {
	this->_AtomList.push_back(lAtoms.getAtom());
    }
    this->_Frames.clear();
}


#ifdef XML_ARCHIVE
    void	Trajectory_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("matter",this->_Matter);
    node->attribute("namespace",this->_Namespace);
    node->archiveVector0("atomList",this->_AtomList);
    node->archiveVector0("frames",this->_Frames);
}
#endif



CL_LISPIFY_NAME("addFrame");
CL_DEFMETHOD TrajectoryFrame_sp Trajectory_O::addFrame(Matter_sp matter)
{_G();
    ASSERTP(matter == this->_Matter,"The matter argument must match the Matter used to define this trajectory");
    TrajectoryFrame_sp frame = TrajectoryFrame_O::create();
    frame->fillFromMatter(this->_AtomList);
    this->_Frames.push_back(frame);
    return frame;
}


CL_LISPIFY_NAME("numberOfTrajectoryFrames");
CL_DEFMETHOD uint Trajectory_O::numberOfTrajectoryFrames()
{_G();
    return this->_Frames.size();
}


CL_LISPIFY_NAME("getTrajectoryFrame");
CL_DEFMETHOD TrajectoryFrame_sp Trajectory_O::getTrajectoryFrame(uint i)
{_G();
    ASSERT_lt(i,this->_Frames.size());
    return this->_Frames[i];
}

#if 0
geom::Render_sp Trajectory_O::rendered(core::List_sp kargs)
{_G();
    geom::FrameList_sp frames = geom::FrameList_O::create();
    for ( gctools::Vec0<TrajectoryFrame_sp>::iterator it=this->_Frames.begin(); it!=this->_Frames.end(); it++ )
    {
	this->applyTrajectoryFrameToMatter(*it);
	frames->append(this->_Matter->rendered(kargs));
    }
    return frames;
}
#endif

CL_LISPIFY_NAME("applyTrajectoryFrameToMatter");
CL_DEFMETHOD void Trajectory_O::applyTrajectoryFrameToMatter(TrajectoryFrame_sp f)
{_G();
    f->applyToMatter(this->_AtomList);
}



};



