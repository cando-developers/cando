/*
    File: dofType.cc
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


#include <clasp/core/foundation.h>
#include <cando/kinematics/dofType.h>

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
