#ifndef	openmmExpose_H //[
#define openmmExpose_H

#include "openmm/include/OpenMM.h"


#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/foundation.h"
#include "core/object.h"
#include "core/ovector3.h"
#include "core/vector3.h"
#include "core/externalObject.h"
#include "core/coordinateArray.h"
#include "units/quantity.fwd.h"
#include "chem/atom.fwd.h"
#include "chem/bond.fwd.h"
#include "chem/angle.fwd.h"
#include "chem/properTorsion.fwd.h"
#include "chem/improperTorsion.fwd.h"
#include "chem/forceField.fwd.h"



namespace omm
{

    static const string ParticleIndex = "ParticleIndex";



#define omm_DeclarePredefinedSymbols
#include "openmmExpose.scrape.inc"




    SMART(System);
    SMART(Integrator);
    SMART(Platform);
    SMART(State);

    template <class ommT>
	core::MetaClass_sp metaClassForOmmClass(core::Lisp_sp lisp)
    {_G();
	THROW(lisp->error(BF("You must specialize template function metaClassForOmmClass (in openmmExpose.cc) for the OpenMM class[%s]") % typeid(ommT).name()  ));
    }



    core::ExternalObject_sp wrap_ommObject(void* ptr, core::MetaClass_sp mc, core::Lisp_sp lisp);




#if 0

template <>
struct	from_object<wxString>
{
    typedef	wxString	ExpectedType;
    typedef	wxString	DeclareType;
    static ExpectedType convert(core::T_sp o)
    {_G();
	if ( o->isAssignableTo<Text_O>() )
	{
	    Text_sp oi = o->as<Text_O>();
	    return wxString(oi->get().c_str());
	}
	THROW(LispError_O::create(BF("Cannot convert %s to wxString") % o->className() ,_lisp));
    }
};

template <>
struct	to_object<const wxString&>
{
    typedef	core::Str_sp		ExpectedType;
    static core::T_sp convert(core::Lisp_sp e,const wxString& v)
    {
	core::Str_sp oi = core::Str_O::create(v.c_str(),e);
	return oi;
    }
};


template <>
struct	to_object<wxString>
{
    typedef	core::Str_sp		ExpectedType;
    static core::T_sp convert(core::Lisp_sp e,const wxString& v)
    {
	core::Str_sp oi = core::Str_O::create(v.c_str(),e);
	return oi;
    }
};

#endif









#define	omm_header
#include "openmmExpose.scrape.inc"
#undef omm_header

};


namespace omm
{

    extern int wrap_Platform_getNumPlatforms(core::Lisp_sp lisp);
    extern Platform_sp wrap_Platform_getPlatform(int index, core::Lisp_sp lisp);
    extern Platform_sp wrap_Platform_getPlatformByName(string const& name, core::Lisp_sp lisp);
    extern core::Cons_sp wrap_Platform_loadPluginsFromDirectory(string const& directoryName, core::Lisp_sp lisp);
};




namespace translate
{
    template <>
    struct	to_object<OpenMM::Vec3>
    {
	typedef	core::OVector3_sp	ExpectedType;
	static core::T_sp convert(core::Lisp_sp l, const OpenMM::Vec3& v)
	{_G();
	    core::OVector3_sp ov = core::OVector3_O::create(l);
	    ov->setAll(v[0],v[1],v[2]);
	    return ov;
	}
    };

    template <>
    struct	to_object<const vector< OpenMM::Vec3>& >
    {
	typedef	core::CoordinateArray_sp	ExpectedType;
	static core::CoordinateArray_sp convert(core::Lisp_sp l, const vector<OpenMM::Vec3>& v)
	{_G();
	    uint numVectors = v.size();
	    core::CoordinateArray_sp ca = core::CoordinateArray_O::create(numVectors,_lisp);
	    for ( uint ii=0; ii<numVectors; ii++ )
	    {
		Vector3& cav = ca->getElement(ii);
		cav.set(v[ii][0],v[ii][1],v[ii][2]);
	    }
	    return ca;
	}
    };





    template <>
    struct	from_object<const OpenMM::Vec3&>
    {
	typedef	OpenMM::Vec3	ExpectedType;
	typedef	OpenMM::Vec3	DeclareType;
	static ExpectedType convert(core::T_sp o)
	{_G();
	    if ( o->isAssignableTo<core::OVector3_O>() )
	    {
		core::OVector3_sp ov = core::OVector3_O::create(_lisp);
		ExpectedType vv;
		vv[0] = ov->getX();
		vv[1] = ov->getY();
		vv[2] = ov->getZ();
		return vv;
	    }
	    THROW(core::LispError_O::create(BF("Cannot convert %s to OpenMM:Vec3") % o->className() ,_lisp));
	}
    };



    template <>
    struct	from_object<const vector<OpenMM::Vec3 >&>
    {
	typedef	vector<OpenMM::Vec3>	ExpectedType;
	typedef	vector<OpenMM::Vec3>	DeclareType;
	static ExpectedType convert(core::T_sp o)
	{_G();
	    if ( o->isAssignableTo<core::CoordinateArray_O>() )
	    {
		core::CoordinateArray_sp ca = o->as<core::CoordinateArray_O>();
		DeclareType res(ca->size());
		for ( int ii=0; ii<(int)ca->size(); ii++ )
		{
		    Vector3& v = ca->getElement(ii);
		    res[ii][0] = v.getX();
		    res[ii][1] = v.getY();
		    res[ii][2] = v.getZ();
		}
		return res;
	    }
	    THROW(core::LispError_O::create(BF("Cannot convert %s to vector<OpenMM:Vec3>") % o->className() ,_lisp));
	}
    };


};
#define	omm_Translate
#include "openmmExpose.scrape.inc"
#undef omm_Translate



#endif //]
