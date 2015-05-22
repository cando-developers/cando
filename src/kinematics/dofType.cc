

#include "core/foundation.h"
#include "dofType.h"

namespace kinematics
{
    string DofType::asString() const
    {
	switch (_Dof)
	{
	case distance: return "distance";
	case theta:	return "theta";
	case phi:	return "phi";
	case rbx:	return "rbx";
	case rby:	return "rby";
	case rbz:	return "rbz";
	case rbalpha:	return "rbalpha";
	case rbbeta:	return "rbbeta";
	case rbgamma:	return "rbgamma";
	default:
	    return "dof-unknown";
	};
    }


    bool DofType::isRigidBodyDof() const
    {
	switch (_Dof)
	{
	case distance: 
	case theta:	
	case phi:	
	    return false;
	case rbx:
	case rby:
	case rbz:
	case rbalpha:
	case rbbeta:
	case rbgamma:
	    return true;
	default:
	    return false;
	};
    }

	    

	    
};

