/*
    File: boundingBox.cc
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
#include <cando/geom/boundingBox.h>
#include <clasp/core/lisp.h>
#include <clasp/core/wrappers.h>



namespace geom 
{


    Vector3 BoundingBox_O::getCenter() const
    {
	Vector3 acc = this->_MinCorner;
	acc = acc + this->_MaxCorner;
	acc = acc * 0.5;
	return acc;
    }

    double BoundingBox_O::getExtentX() const { return this->_MaxCorner.getX()-this->_MinCorner.getX();};
    double BoundingBox_O::getExtentY() const { return this->_MaxCorner.getY()-this->_MinCorner.getY();};
    double BoundingBox_O::getExtentZ() const { return this->_MaxCorner.getZ()-this->_MinCorner.getZ();};



CL_LISPIFY_NAME("getMinX");
CL_DEFMETHOD     double BoundingBox_O::getMinX() const { return this->_MinCorner.getX(); };
CL_LISPIFY_NAME("getMaxX");
CL_DEFMETHOD     double BoundingBox_O::getMaxX() const { return this->_MaxCorner.getX(); };
CL_LISPIFY_NAME("getMinY");
CL_DEFMETHOD     double BoundingBox_O::getMinY() const { return this->_MinCorner.getY(); };
CL_LISPIFY_NAME("getMaxY");
CL_DEFMETHOD     double BoundingBox_O::getMaxY() const { return this->_MaxCorner.getY(); };
CL_LISPIFY_NAME("getMinZ");
CL_DEFMETHOD     double BoundingBox_O::getMinZ() const { return this->_MinCorner.getZ(); };
CL_LISPIFY_NAME("getMaxZ");
CL_DEFMETHOD     double BoundingBox_O::getMaxZ() const { return this->_MaxCorner.getZ(); };



    void	BoundingBox_O::initialize()
    {
	this->Base::initialize();
	this->_Defined = false;
	this->_MinCorner.set(0.0,0.0,0.0);
	this->_MaxCorner.set(0.0,0.0,0.0);
    }

#if defined(XML_ARCHIVE)
    void	BoundingBox_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attributeIfNotDefault("defined",this->_Defined,true);
	node->archivePlainObject<Vector3>("minCorner","MinCorner",this->_MinCorner);
	node->archivePlainObject<Vector3>("maxCorner","MaxCorner",this->_MaxCorner);
    }
#endif

    void	BoundingBox_O::expandToEncompassPoint(const Vector3& p)
    {
	if ( !this->_Defined )
	{
	    this->_Defined = true;
	    this->_MinCorner = p;
	    this->_MaxCorner = p;
	    return;
	}
	double minX = MIN(p.getX(),this->_MinCorner.getX());
	double minY = MIN(p.getY(),this->_MinCorner.getY());
	double minZ = MIN(p.getZ(),this->_MinCorner.getZ());
	double maxX = MAX(p.getX(),this->_MaxCorner.getX());
	double maxY = MAX(p.getY(),this->_MaxCorner.getY());
	double maxZ = MAX(p.getZ(),this->_MaxCorner.getZ());
	this->_MinCorner.set(minX,minY,minZ);
	this->_MaxCorner.set(maxX,maxY,maxZ);
    }


    void BoundingBox_O::expandToEncompassPoint(const Matrix& transform, const Vector3& point )
    {
	Vector3 tpnt = transform*point;
	this->expandToEncompassPoint(tpnt);
    }


    bool	BoundingBox_O::encompasses(const Vector3& p)
    {_OF();
	ASSERT(this->_Defined);
	if ( !(this->_MinCorner.getX() <= p.getX() && p.getX() <= this->_MaxCorner.getX() ) ) return false;
	if ( !(this->_MinCorner.getY() <= p.getY() && p.getY() <= this->_MaxCorner.getY() ) ) return false;
	if ( !(this->_MinCorner.getZ() <= p.getZ() && p.getZ() <= this->_MaxCorner.getZ() ) ) return false;
	return true;
    }


CL_LISPIFY_NAME("pad");
CL_DEFMETHOD     void BoundingBox_O::pad(double add)
    {_OF();
	ASSERT(this->_Defined);
	Vector3 vadd;
	vadd.set(add,add,add);
	this->_MinCorner = this->_MinCorner.sub(vadd);
	this->_MaxCorner = this->_MaxCorner.add(vadd);
    }






};
