/*
    File: calculatePosition.cc
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
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/residue.h>
#include <clasp/core/wrappers.h>


namespace chem {



void	CalculatePosition_O::initialize()
{
    this->Base::initialize();
}

#ifdef XML_ARCHIVE
    void	CalculatePosition_O::archiveBase(core::ArchiveP node)
{
    // nothing
}
#endif

#if INIT_TO_FACTORIES

#define ARGS_CalculatePositionUsingInternals_O_make "(distance_atom_name angle_atom_name dihedral_atom_name distance angle_degrees dihedral_degrees)"
#define DECL_CalculatePositionUsingInternals_O_make ""
#define DOCS_CalculatePositionUsingInternals_O_make "make AtomId"
  CalculatePositionUsingInternals_sp CalculatePositionUsingInternals_O::make(MatterName distanceAtomName, MatterName angleAtomName, MatterName dihedralAtomName,
									     double distance, double angleDegrees, double dihedralDegrees)
  {
      GC_ALLOCATE(CalculatePositionUsingInternals_O, me );
    me->_DistanceAtomName = distanceAtomName;
    me->_AngleAtomName = angleAtomName;
    me->_DihedralAtomName = dihedralAtomName;
    me->_Distance = distance;
    me->_AngleDegrees = angleDegrees;
    me->_DihedralDegrees = dihedralDegrees;
    return me;
  };

#else

    core::T_sp 	CalculatePositionUsingInternals_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{
    this->Base::__init__(exec,args,env,lisp);
    this->_DistanceAtomName = env->lookup(lisp->internWithPackageName(ChemPkg,"distanceAtomName")).as<core::Str_O>()->get();
    this->_AngleAtomName = env->lookup(lisp->internWithPackageName(ChemPkg,"angleAtomName")).as<core::Str_O>()->get();
    this->_DihedralAtomName = env->lookup(lisp->internWithPackageName(ChemPkg,"dihedralAtomName")).as<core::Str_O>()->get();
    this->_Distance = env->lookup(lisp->internWithPackageName(ChemPkg,"distance")).as<core::DoubleFloat_O>()->get();
    this->_AngleDegrees = env->lookup(lisp->internWithPackageName(ChemPkg,"angleDegrees")).as<core::DoubleFloat_O>()->get();
    this->_DihedralDegrees = env->lookup(lisp->internWithPackageName(ChemPkg,"dihedralDegrees")).as<core::DoubleFloat_O>()->get();
    return _Nil<core::T_O>();
}

#endif

void	CalculatePositionUsingInternals_O::initialize()
{
    this->Base::initialize();
    this->_DistanceAtomName = _Nil<MatterName::Type>();
    this->_AngleAtomName = _Nil<MatterName::Type>();
    this->_DihedralAtomName = _Nil<MatterName::Type>();
}

#ifdef XML_ARCHIVE
    void	CalculatePositionUsingInternals_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("DistanceAtomName",this->_DistanceAtomName);
    node->attribute("AngleAtomName",this->_AngleAtomName);
    node->attribute("DihedralAtomName",this->_DihedralAtomName);
    node->attribute("Distance",this->_Distance);
    node->attribute("AngleDegrees",this->_AngleDegrees);
    node->attribute("DihedralDegrees",this->_DihedralDegrees);
}
#endif

    Vector3 CalculatePositionUsingInternals_O::calculatePosition(core::T_sp positionLookup)
{
	// Calculate the position according to the arguments
	//
    Vector3 pnt,posDistanceAtom,posAngleAtom,posDihedralAtom;

#if !KINEMATICS
    if ( positionLookup.isA<IterateFragmentCoordinateEntries_O>() )
    {
	IterateFragmentCoordinateEntries_sp it = positionLookup.as<IterateFragmentCoordinateEntries_O>();
	posDistanceAtom = it->positionOfAtomWithName(this->_DistanceAtomName);
	posAngleAtom = it->positionOfAtomWithName(this->_AngleAtomName);
	posDihedralAtom = it->positionOfAtomWithName(this->_DihedralAtomName);
    } else
#endif
	if ( positionLookup.isA<Residue_O>() )
    {
	Residue_sp res = positionLookup.as<Residue_O>();
	posDistanceAtom = res->positionOfAtomWithName(this->_DistanceAtomName);
	posAngleAtom = res->positionOfAtomWithName(this->_AngleAtomName);
	posDihedralAtom = res->positionOfAtomWithName(this->_DihedralAtomName);
    } else
    {
	SIMPLE_ERROR(BF("Unacceptable positionLookup class"));
    }
    LOG(BF("About to calculate point with values") );
    LOG(BF(" distance(%lf)  distancePoint(%s)") % this->_Distance % posDistanceAtom.asString().c_str()  );
    LOG(BF(" angleDegrees(%lf)  anglePoint(%s)") % this->_AngleDegrees % posAngleAtom.asString().c_str()  );
    LOG(BF(" dihedralDegrees(%lf)  dihedralPoint(%s)") % this->_DihedralDegrees % posDihedralAtom.asString().c_str()  );
    pnt = geom::buildUsingBondAngleDihedral(this->_Distance,posDistanceAtom,
			this->_AngleDegrees*0.0174533,posAngleAtom,
				      this->_DihedralDegrees*0.0174533,posDihedralAtom);
    LOG(BF("Generated point: %s") % pnt.asString().c_str()  );
    return pnt;
}

#if INIT_TO_FACTORIES

#define ARGS_CalculatePositionAlongBond_O_make "(from_atom_name to_atom_name extend_distance)"
#define DECL_CalculatePositionAlongBond_O_make ""
#define DOCS_CalculatePositionAlongBond_O_make "make CalculatePositionAlongBond"
  CalculatePositionAlongBond_sp CalculatePositionAlongBond_O::make(MatterName fromAtomName, MatterName toAtomName, double extendDistance)
  {
      GC_ALLOCATE(CalculatePositionAlongBond_O, me );
    me->_FromAtomName = fromAtomName;
    me->_ToAtomName = toAtomName;
    me->_ExtendDistance = extendDistance;
    return me;
  };

#else

    core::T_sp 	CalculatePositionAlongBond_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{
    this->Base::__init__(exec,args,env,lisp);
    this->_FromAtomName = env->lookup(lisp->internWithPackageName(ChemPkg,"fromAtomName")).as<core::Str_O>()->get();
    this->_ToAtomName = env->lookup(lisp->internWithPackageName(ChemPkg,"toAtomName")).as<core::Str_O>()->get();
    this->_ExtendDistance = env->lookup(lisp->internWithPackageName(ChemPkg,"extendDistance")).as<core::DoubleFloat_O>()->get();
    return _Nil<core::T_O>();
}

#endif

void	CalculatePositionAlongBond_O::initialize()
{
    this->Base::initialize();
    this->_FromAtomName = _Nil<MatterName::Type>();
    this->_ToAtomName = _Nil<MatterName::Type>();
    this->_ExtendDistance = 0.0;
}

#ifdef XML_ARCHIVE
    void	CalculatePositionAlongBond_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("FromAtomName",this->_FromAtomName);
    node->attribute("ToAtomName",this->_ToAtomName);
    node->attribute("ExtendDistance",this->_ExtendDistance);
}
#endif




    Vector3 CalculatePositionAlongBond_O::calculatePosition(core::T_sp positionLookup)
{
    IMPLEMENT_ME();
#if !KINEMATICS
	// Calculate the position according to the arguments
	//
    Vector3 pnt,posFromAtom,posToAtom;
    	// If the 

    if ( positionLookup.isA<IterateFragmentCoordinateEntries_O>() )
    {
	IterateFragmentCoordinateEntries_sp it = positionLookup.as<IterateFragmentCoordinateEntries_O>();
	posFromAtom = it->positionOfAtomWithName(this->_FromAtomName);
	posToAtom = it->positionOfAtomWithName(this->_ToAtomName);
	LOG(BF("Got atom positions from IterateFragmentCoordinateEntries with names (%s)-(%s)") % this->_FromAtomName.c_str() % this->_ToAtomName.c_str()  );
    } else if ( positionLookup.isA<Residue_O>() )
    {
	Residue_sp res = positionLookup.as<Residue_O>();
	posFromAtom = res->positionOfAtomWithName(this->_FromAtomName);
	posToAtom = res->positionOfAtomWithName(this->_ToAtomName);
	LOG(BF("Got atom positions from Residue with names (%s)-(%s)") % this->_FromAtomName.c_str() % this->_ToAtomName.c_str()  );
    } else
    {
	SIMPLE_ERROR(BF("Unacceptable positionLookup class"));
    }
    LOG(BF("About to calculate point with values") );
    LOG(BF(" fromAtom(%s)   -  toAtom(%s)") % posFromAtom.asString().c_str() % posToAtom.asString().c_str()  );
    LOG(BF(" extendDistance(%lf) ") % this->_ExtendDistance  );
    Vector3 dir = posToAtom.sub(posFromAtom);
    dir = dir.normalizedOrZero();
    if ( !dir.isDefined() )
    {
		// If the atom positions have not yet
		// been defined then just return 0,0,0
	return dir;
    }
    pnt = posToAtom.add(dir.multiplyByScalar(this->_ExtendDistance));
    LOG(BF("Generated point: %s") % pnt.asString().c_str()  );
    return pnt;
#endif
}

#if INIT_TO_FACTORIES

#define ARGS_CalculatePositionRelativeToOrigin_O_make "(origin_atom_name to_atom_name extend_distance)"
#define DECL_CalculatePositionRelativeToOrigin_O_make ""
#define DOCS_CalculatePositionRelativeToOrigin_O_make "make CalculatePositionRelativeToOrigin"
  CalculatePositionRelativeToOrigin_sp CalculatePositionRelativeToOrigin_O::make(MatterName originAtomName, MatterName toAtomName, double extendDistance)
  {
      GC_ALLOCATE(CalculatePositionRelativeToOrigin_O, me );
    me->_OriginAtomName = originAtomName;
    me->_ToAtomName = toAtomName;
    me->_ExtendDistance = extendDistance;
    return me;
  };

#else

    core::T_sp 	CalculatePositionRelativeToOrigin_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{
    this->Base::__init__(exec,args,env,lisp);
    this->_OriginAtomName = env->lookup(lisp->internWithPackageName(ChemPkg,"originAtomName")).as<core::Str_O>()->get();
    this->_ToAtomName = env->lookup(lisp->internWithPackageName(ChemPkg,"toAtomName")).as<core::Str_O>()->get();
    this->_ExtendDistance = env->lookup(lisp->internWithPackageName(ChemPkg,"extendDistance")).as<core::DoubleFloat_O>()->get();
    return _Nil<core::T_O>();
}

#endif

void	CalculatePositionRelativeToOrigin_O::initialize()
{
    this->Base::initialize();
    this->_OriginAtomName = _Nil<MatterName::Type>();
    this->_ToAtomName = _Nil<MatterName::Type>();
    this->_ExtendDistance = 0.0;
}

#ifdef XML_ARCHIVE
    void	CalculatePositionRelativeToOrigin_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("OriginAtomName",this->_OriginAtomName);
    node->attribute("ToAtomName",this->_ToAtomName);
    node->attribute("ExtendDistance",this->_ExtendDistance);
}
#endif



/*
 * Sometimes the fragment won't contain one of the atoms that we want to
 * use to calculate a VirtualAtom but that atom will be on the origin so
 * we use that.
 */

    Vector3 CalculatePositionRelativeToOrigin_O::calculatePosition(core::T_sp positionLookup)
{
    IMPLEMENT_ME();
#if !KINEMATICS
	// Calculate the position according to the arguments
	//
    Vector3 pnt,posOriginAtom,posToAtom;
    	// If the 
    if ( positionLookup.isA<IterateFragmentCoordinateEntries_O>() )
    {
	IterateFragmentCoordinateEntries_sp it = positionLookup.as<IterateFragmentCoordinateEntries_O>();
	posOriginAtom.set(0.0,0.0,0.0);
	posToAtom = it->positionOfAtomWithName(this->_ToAtomName);
	LOG(BF("Got atom positions from IterateFragmentCoordinateEntries with names (%s)-(%s)") % this->_OriginAtomName.c_str() % this->_ToAtomName.c_str()  );
    } else if ( positionLookup.isA<Residue_O>() )
    {
	Residue_sp res = positionLookup.as<Residue_O>();
	posOriginAtom = res->positionOfAtomWithName(this->_OriginAtomName);
	posToAtom = res->positionOfAtomWithName(this->_ToAtomName);
	LOG(BF("Got atom positions from Residue with names (%s)-(%s)") % this->_OriginAtomName.c_str() % this->_ToAtomName.c_str()  );
    } else
    {
	SIMPLE_ERROR(BF("Unacceptable positionLookup class"));
    }
    LOG(BF("About to calculate point with values") );
    LOG(BF(" originAtom(%s)   -  toAtom(%s)") % posOriginAtom.asString().c_str() % posToAtom.asString().c_str()  );
    LOG(BF(" extendDistance(%lf) ") % this->_ExtendDistance  );
    Vector3 dir = posToAtom.sub(posOriginAtom);
    dir = dir.normalizedOrZero();
    if ( !dir.isDefined() )
    {
		// If the atom positions have not yet
		// been defined then just return 0,0,0
	return dir;
    }
    pnt = posToAtom.add(dir.multiplyByScalar(this->_ExtendDistance));
    LOG(BF("Generated point: %s") % pnt.asString().c_str()  );
    return pnt;
#endif
}
































};
