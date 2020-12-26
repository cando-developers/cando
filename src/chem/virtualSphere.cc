/*
    File: virtualSphere.cc
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


//
// (C) 2004 Christian E. Schafmeister
//


#include <stdio.h>
#include <math.h>
#include <clasp/core/common.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/corePackage.h>
#include <cando/chem/virtualSphere.h>
#include <clasp/core/symbolTable.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>


#define	SCALE_IN_VIRTUAL_SPHERE 1


namespace chem 
{

    SYMBOL_EXPORT_SC_(ChemPkg,idle);
    SYMBOL_EXPORT_SC_(ChemPkg,translate);
    SYMBOL_EXPORT_SC_(ChemPkg,rotate);
    SYMBOL_EXPORT_SC_(ChemPkg,scale);


/*
 * matrix_ident:
 *	Create an identity matrix
 */
static void matrix_identity( double m[4][4] )
{
  int	i, j;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++)
			m[i][j] = 0;
		m[i][i] = 1;
	}
}

/*
 * matrix_multiply:
 *	Multiply two matrices
 */
static void
matrix_multiply(double m1[4][4], double m2[4][4], Matrix& result )
{
   int	i, j, k;
	double		tmp;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			tmp = 0;
			for (k = 0; k < 4; k++)
				tmp += m1[i][k] * m2[k][j];
			result.atRowColPut(i,j,tmp);
		}
}

/*
 * vxv3:
 *	Take the cross product of two vectors
 */
static void
vxv3(double c[3], double a[3], double b[3])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

/*
 * normalize:
 *	Normalize a vector
 */
static void
normalize(double v[3])
{
double	d;

	d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v[0] /= d;
	v[1] /= d;
	v[2] /= d;
}


/*
 * vsphere:
 *	Convert X-Y coordinates into a 4x4 rotation matrix
 */
int vsphere(double fx, double fy, double tx, double ty, Matrix& rotmat )
{
   int	i;
	double		d1, d2;
	double		f[3], t[3], a[3], g[3], u[3];
	double		m1[4][4], m2[4][4];

	/*
	 * First construct the unit vectors and computing
	 * the normal unit vectors.  If we cross from outside
	 * the sphere to inside the sphere (or vice versa), then
	 * we ignore the transition because the discontinuity
	 * tends to make molecules jump.
	 */
	d1 = fx * fx + fy * fy;
	d2 = tx * tx + ty * ty;
	if (d1 > 1 && d2 < 1)
		return -2;
	if (d2 > 1 && d1 < 1)
		return -3;
	if (d1 < 1) {
		f[0] = fx;
		f[1] = fy;
		f[2] = sqrt(1 - d1);
	}
	else {
		d1 = sqrt(d1);
		f[0] = fx / d1;
		f[1] = fy / d1;
		f[2] = 0;
	}
	if (d2 < 1) {
		t[0] = tx;
		t[1] = ty;
		t[2] = sqrt(1 - d2);
	}
	else {
		d2 = sqrt(d2);
		t[0] = tx / d2;
		t[1] = ty / d2;
		t[2] = 0;
	}

	/*
	 * If the positions normalize to the same place we just punt.
	 * We don't even bother to put in the identity matrix.
	 */
	if (f[0] == t[0] && f[1] == t[1] && f[2] == t[2])
		return -1;

	vxv3(a, f, t);
	normalize(a);
	vxv3(g, a, f);		/* Don't need to normalize these since the */
	vxv3(u, a, t);		/* cross product of normal unit vectors is */
				/* a unit vector */

	/*
	 * Now assemble them into the inverse matrix (to go from
	 * the from-vector to xyz-space) and the transform matrix
	 * (to go from xyz-space to to-vector).  The product of
	 * the inverse and transformation matrices is the rotation
	 * matrix.
	 */
	matrix_identity(m1);
	matrix_identity(m2);
	for (i = 0; i < 3; i++) {
		m1[i][0] = f[i];
		m1[i][1] = a[i];
		m1[i][2] = g[i];
		m2[0][i] = t[i];
		m2[1][i] = a[i];
		m2[2][i] = u[i];
	}
	matrix_multiply(m1, m2, rotmat);
	return 0;
}


// ##############################################################
// ##############################################################
// ##############################################################
// ##############################################################
// ##############################################################











void	VirtualSphere_O::convertToRelative( int ax, int ay,
					double& rx, double& ry )
{
int	halfWidth, halfHeight;
double	radius;
    halfWidth = this->_WindowWidth/2;
    halfHeight = this->_WindowHeight/2;
    radius = fmin(halfWidth,halfHeight)*0.8;
    rx = (ax-halfWidth)/radius;
//    ry = (halfHeight - ay)/radius;
//    rx = (halfWidth-ax)/radius;
    ry = (ay-halfHeight)/radius;
//    LOG(BF("relative mouse position = %lf, %lf") % (rx) % (ry ) );
}







#ifdef XML_ARCHIVE
    void	VirtualSphere_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("windowHeight",this->_WindowHeight);
    node->attribute("windowWidth",this->_WindowWidth);
    node->attribute("startX",this->_StartX);
    node->attribute("startY",this->_StartY);
    node->attribute("mouseX",this->_MouseX);
    node->attribute("mouseY",this->_MouseY);
//    node->attribute("mouseMode",this->_MouseMode);
    node->attribute("recalculate",this->_Recalculate);
    node->attribute("scale",this->_Scale);
    node->archivePlainObject<Vector3>( "translate","Vector3",this->_Translate);
    node->archivePlainObject<Vector3>( "center","Vector3",this->_Center);
    node->archivePlainObject( "rotate","Matrix",this->_Rotate);
    node->archivePlainObject( "overall","Matrix",this->_Overall);
}
#endif


void	VirtualSphere_O::initialize()
{
    this->Base::initialize();
    this->_MouseMode = vsIdle;
    this->_Recalculate = true;
    this->_Scale = 0.25;
    this->_Translate.set(0.0,0.0,0.0);
    this->_Center.set(0.0,0.0,0.0);
}


CL_LISPIFY_NAME("dump");
CL_DEFMETHOD void VirtualSphere_O::dump()
{_OF();
    stringstream ss;
    ss << "VirtualSphere dump  _Recalculate = " << this->_Recalculate << std::endl;
    ss << "_Scale = " << this->_Scale << std::endl;
    ss << "_Translate = " << this->_Translate.asString() << std::endl;
    ss << "_Center = " << this->_Center.asString() << std::endl;
    ss << "_Rotate = " << this->_Rotate.asString() << std::endl;
    printf("%s", ss.str().c_str() );
    LOG(BF("%s") % ss.str() );
}

SYMBOL_EXPORT_SC_(ChemKwPkg,idle);
SYMBOL_EXPORT_SC_(ChemKwPkg,translate);
SYMBOL_EXPORT_SC_(ChemKwPkg,rotate);
SYMBOL_EXPORT_SC_(ChemKwPkg,scale);

    VirtualSphereMotionMode 
    VirtualSphere_O::mouseModeFromSymbol(core::Symbol_sp sym)
    {_OF();
	VirtualSphereMotionMode mode = vsIdle;
	if ( sym == chemkw::_sym_idle ) {mode = vsIdle;}
	else if ( sym == chemkw::_sym_translate ) {mode = vsTranslate;}
	else if ( sym == chemkw::_sym_rotate ) {mode = vsRotate;}
	else if ( sym == chemkw::_sym_scale ) {mode = vsScale;}
	return mode;
    }


CL_LISPIFY_NAME("setMouseStart");
CL_DEFMETHOD     void VirtualSphere_O::setMouseStartSymbol( core::Symbol_sp sym, int x, int y, int w, int h )
    {_OF();
	VirtualSphereMotionMode mode = this->mouseModeFromSymbol(sym);
	this->setMouseStart(mode,x,y,w,h);
    }
    void VirtualSphere_O::setMouseStart( VirtualSphereMotionMode  mode, int x, int y, int w, int h )
    {
	this->_MouseMode = mode;
	this->_StartX = x;
	this->_StartY = y;
	this->_MouseX = x;
	this->_MouseY = y;
	this->_WindowHeight = h;
	this->_WindowWidth = w;
	this->_Recalculate = true;
    };



CL_LISPIFY_NAME("setMouseMove");
CL_DEFMETHOD void	VirtualSphere_O::setMouseMove( int x, int y )
{_OF();
double	sx,sy,tx,ty,scaleDelta;
Vector3	delta;
int	result;
Matrix	mtemp;
    this->_MouseX = x;
    this->_MouseY = y;
    switch ( this->_MouseMode )
    {
    case vsIdle:
	break;
    case vsRotate:
	this->convertToRelative( this->_StartX, this->_StartY, sx,sy );
	this->convertToRelative( this->_MouseX, this->_MouseY, tx,ty );
//	    result = vsphere( sx,sy,tx,ty, mtemp.elements );
	result = vsphere( tx,ty,sx,sy, mtemp );
	if ( result == 0 ) {
	    this->_Rotate = mtemp.multiplyByMatrix(this->_Rotate);
//	        this->_Rotate = this->_Rotate.multiplyByMatrix(mtemp);
	}
	this->_StartX = this->_MouseX;
	this->_StartY = this->_MouseY;
	this->_Recalculate = true;
        FIX_ME(); //this->notify(chemkw::_sym_changed,_Nil<core::T_O>());
	break;
    case vsTranslate:
	this->convertToRelative( this->_StartX, this->_StartY, sx,sy );
	this->convertToRelative( this->_MouseX, this->_MouseY, tx,ty );
	delta.set(tx-sx,ty-sy,0.0);
	this->_Translate = this->_Translate.add(delta*10);
//	    printf("translation delta = %f,%f,%f\n", tx-sx,ty-sy,0.0 );
//	    printf("translation = %f,%f,%f\n",
//				this->_Translate.getX(),
//				this->_Translate.getY(),
//				this->_Translate.getZ() );
	this->_StartX = this->_MouseX;
	this->_StartY = this->_MouseY;
	this->_Recalculate = true;
	FIX_ME(); // this->notify(chemkw::_sym_changed,_Nil<core::T_O>());
	break;
    case vsScale:
	this->convertToRelative( this->_StartX, this->_StartY, sx,sy );
	this->convertToRelative( this->_MouseX, this->_MouseY, tx,ty );
	scaleDelta = sy - ty;

//            cout<<"The scale delta is "<<scaleDelta<<"\t";
//	    printf( "Scale delta = %f\n", scaleDelta );
	scaleDelta = pow(2.0,scaleDelta);
//	    printf( "2^(Scale delta) = %f\n", scaleDelta );
	this->_Scale *= scaleDelta;
//	    printf( "this->_Scale = %f\n", this->_Scale);
//	    this->_Scale = MAX(this->_Scale,0.005);
//	    this->_Scale = MIN(this->_Scale,50.0);
	this->_StartX = this->_MouseX;
	this->_StartY = this->_MouseY;
	this->_Recalculate = true;
	break;
    }
}


    void VirtualSphere_O::setScaleChange(double yy)
    {
	double scaleDelta = yy;
//            cout<<"The scale delta is "<<scaleDelta<<"\t";
//	    printf( "Scale delta = %f\n", scaleDelta );
	scaleDelta = pow(2.0,scaleDelta);
//	    printf( "2^(Scale delta) = %f\n", scaleDelta );
	this->_Scale *= scaleDelta;
//	    printf( "this->_Scale = %f\n", this->_Scale);
//	this->_Scale = MAX(this->_Scale,0.005);
//	this->_Scale = MIN(this->_Scale,50.0);
	this->_Recalculate = true;
    }




CL_LISPIFY_NAME("setMouseXAxisWheelChange");
CL_DEFMETHOD     void
    VirtualSphere_O::setMouseXAxisWheelChange(core::Symbol_sp smode, double w)
    {_OF();
	VirtualSphereMotionMode mode = this->mouseModeFromSymbol(smode);
	LOG(BF("Mouse mode=%d") % mode );
	switch ( mode )
	{
	case vsIdle:
	    LOG(BF("Mouse mode=== idle"));
	    break;
	case vsRotate:
	{
	    LOG(BF("Mouse mode=== rotate"));
	    Matrix mtemp;
	    if ( w > 1.0 ) w = 1.0;
	    if ( w < -1.0 ) w = -1.0;
	    mtemp.rotationY(w);
	    this->_Rotate = mtemp.multiplyByMatrix(this->_Rotate);
	    LOG(BF("Rotating with matrix: %s") % (mtemp.asString()));
//	        this->_Rotate = this->_Rotate.multiplyByMatrix(mtemp);
	    this->_Recalculate = true;
	    break;
	}
	case vsTranslate:
	{
	    LOG(BF("Mouse mode=== translate"));
	    Vector3 delta;
	    delta.set(w,0.0,0.0);
	    this->_Translate = this->_Translate.add(delta);
	    this->_Recalculate = true;
	    break;
	}
	case vsScale:
	    LOG(BF("Mouse mode=== scale"));
	    this->setScaleChange(w);
	    break;
    default:
	LOG(BF("Unknown mouse mode"));
	break;
	}
    }


CL_LISPIFY_NAME("setMouseYAxisWheelChange");
CL_DEFMETHOD     void
    VirtualSphere_O::setMouseYAxisWheelChange(core::Symbol_sp smode, double w)
    {_OF();
	VirtualSphereMotionMode mode = this->mouseModeFromSymbol(smode);
	LOG(BF("Mouse mode=%d") % mode );
	switch ( mode )
	{
	case vsIdle:
	    break;
	case vsRotate:
	{
	    LOG(BF("Mouse mode=== rotate"));
	    Matrix mtemp;
	    if ( w > 1.0 ) w = 1.0;
	    if ( w < -1.0 ) w = -1.0;
	    mtemp.rotationX(w);
	    this->_Rotate = mtemp.multiplyByMatrix(this->_Rotate);
	    LOG(BF("Rotating with matrix: %s") % (mtemp.asString()));
//	        this->_Rotate = this->_Rotate.multiplyByMatrix(mtemp);
	    this->_Recalculate = true;
	    break;
	}
	case vsTranslate:
	{
	    LOG(BF("Mouse mode=== translate"));
	    Vector3 delta;
	    delta.set(0.0,w,0.0);
	    this->_Translate = this->_Translate.add(delta);
	    this->_Recalculate = true;
	    break;
	}
	case vsScale:
	    LOG(BF("Mouse mode=== scale"));
	    this->setScaleChange(w);
	    break;
	default:
	    LOG(BF("Unknown mouse mode"));
	    break;
	}
    }



void	VirtualSphere_O::setCenter(const Vector3& pos)
{
    LOG(BF("setCenter to %s") % pos.asString() );
    this->_Center = pos;
    this->_Recalculate = true;
}


void	VirtualSphere_O::setTranslate(const Vector3& pos )
{
    LOG(BF("setTranslate to %s") % pos.asString() );
    this->_Translate = pos;
    this->_Recalculate = true; 
}

/*! Get the scale of the VirtualSphere,
 * do not let it get larger than maxScale or smaller than minScale
 */
double	VirtualSphere_O::getScale(double minScale, double maxScale)
{
    if ( this->_Scale < minScale )
    {
	this->_Scale = minScale;
    }
    if ( this->_Scale > maxScale )
    {
	this->_Scale = maxScale;
    }
    return this->_Scale;
}



Matrix	VirtualSphere_O::getOverallMatrix()
{
Matrix	result;
Vector3	vCent;

    if (this->_Recalculate == false ) {
        LOG(BF("Returning pre-calculated matrix") );
	result = this->_Overall;
	return result;
    }
    LOG(BF("Recalculating matrix") );
#if 1 //[
Matrix	mTranslateCenterToOrigin;
Matrix	mTranslateOriginToCenter;
Matrix	mTranslate;
Matrix	mScale;
Matrix	mA, mB;


		//
		// Calculate matrices
		//
    vCent = this->_Center.multiplyByScalar(-1.0);
    LOG(BF("vCent = %s") % vCent.asString().c_str()  );
    mTranslateCenterToOrigin.translate(vCent);
    mTranslateOriginToCenter.translate(this->_Center);
    mTranslate.translate(this->_Translate);

    mA = this->_Rotate.multiplyByMatrix(mTranslateCenterToOrigin);
#if SCALE_IN_VIRTUAL_SPHERE //[ // Don't scale here, do it in glPerspective
    mScale.scale(this->_Scale);
    mB = mScale.multiplyByMatrix(mA);
#else //] [
    mB = mA;
#endif //]
    mA = mTranslateOriginToCenter.multiplyByMatrix(mB);
    result = mTranslate.multiplyByMatrix(mA);
#endif // ]
#if 0	// Try the above
    mA = this->_Rotate.multiplyByMatrix(mTranslateCenterToOrigin);
    mB = mA.multiplyByMatrix(mTranslateOriginToCenter);
    mA = mB.multiplyByMatrix(mScale);
    result = mA.multiplyByMatrix(mTranslate);
#endif
    this->_Overall = result;
    this->_Recalculate = false;
    return result;
}


};
