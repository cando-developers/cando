       
#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//

//
// coordSys.cc
//
#include <clasp/core/common.h>
#include <cando/chem/coordSys.h>
#include <clasp/core/numerics.h>
#include <cando/geom/ovector3.h>
#include <cando/geom/omatrix.h>
//#include "core/archiveNode.h"
#include <cando/geom/color.h>
#include <cando/chem/atom.h>
#include <clasp/core/wrappers.h>

#define	PI		3.1415926535897932
#define	DEGTORAD	0.0174533
#define	PARALLEL_CUTOFF	0.00001


namespace chem {

#if 0
O_GimbalTransform::O_GimbalTransform()
{
}

O_GimbalTransform::~O_GimbalTransform()
{
}




void	O_GimbalTransform::defineForCoordinateSystem(CoordinateSystem_sp coord)
{_G();
    Vector3	labX; labX.set(1.0,0.0,0.0);
    Vector3	labY; labY.set(0.0,1.0,0.0);
    Vector3	labZ; labZ.set(0.0,0.0,1.0);
    Vector3 c2xo = coord->getX();
    Vector3 xycomp = Vector3( c2xo.getX(), c2xo.getY(), 0.0 ).normalized(_lisp);
        LOG(BF("xycomp={%lf,%lf,%lf};")%xycomp.getX()%xycomp.getY()%xycomp.getZ());
	double cosXy = xycomp.dotProduct(labX);
        LOG(BF("cosXy = %lf;")%cosXy);
        if ( xycomp.dotProduct(labY)>0 ) {
	    this->_RotZ = acos(cosXy);
        } else {
	    this->_RotZ = -acos(cosXy);
	}
	LOG(BF("this->_RotZ = %lf degrees") % this->_RotZ/DEGTORAD  );
            //
            // Now undo the effect of rotz on c2.getX()
	    //
        Matrix mrotz;
        mrotz.rightHandedRotationZ(-this->_RotZ);
	Vector3 c2x; c2x = mrotz.multiplyByVector3(c2xo);
	double cosC2x = c2x.dotProduct(labX);
        if ( c2x.dotProduct(labZ) < 0.0 ) {
	    this->_RotY = acos(cosC2x);
        } else {
	    this->_RotY = -acos(cosC2x);
	}
	LOG(BF("this->_RotY = %lf degrees") % this->_RotY/DEGTORAD  );
        Matrix mroty;
        mroty.rightHandedRotationY(-this->_RotY);
	Vector3 c3x = mroty.multiplyByVector3(c2x);
	Vector3 c2y = mrotz.multiplyByVector3(coord->getY());
	Vector3 c3y = mroty.multiplyByVector3(c2y);
    //    renderRay(self._gf,YELLOW,c1.getOrigin(),c3y)
        if ( c3y.dotProduct(labZ) > 0.0 ) {
            this->_RotX = acos(c3y.dotProduct(labY));
        } else {
	    this->_RotX = -acos(c3y.dotProduct(labY));
	}
	LOG(BF("this->_RotX = %lf degrees") % this->_RotX/DEGTORAD  );
	Vector3 trans = coord->getOrigin();
	Vector3 transNorm = trans.normalized(_lisp);
        LOG(BF("trans= {%lf,%lf,%lf};")%trans.getX()%trans.getY()%trans.getZ());
        xycomp = Vector3( trans.getX(), trans.getY(),0.0 ).normalized(_lisp);
        LOG(BF("xycomp= {%lf,%lf,%lf};") % xycomp.getX() % xycomp.getY() % xycomp.getZ() );
	cosXy = xycomp.dotProduct(labX);
	LOG(BF("cosXy = %lf") % cosXy );
        if ( xycomp.dotProduct(labY)>0 )
	{
	    this->_RotZZ = acos(cosXy);
        } else {
	    this->_RotZZ = -acos(cosXy);
	}
	LOG(BF("this->_RotZZ = %lf degrees") % this->_RotZZ/DEGTORAD  );
            //
            // Now undo the effect of rotz on c2.getX()
	    //
        mrotz.rightHandedRotationZ(-this->_RotZZ);
        c2x = mrotz.multiplyByVector3(transNorm);
        LOG(BF("c2x= {%lf,%lf,%lf};") % c2x.getX() % c2x.getY() % c2x.getZ() );
	cosC2x = c2x.dotProduct(labX);
	LOG(BF("cosC2x = %lf") % cosC2x  );
        if ( c2x.dotProduct(labZ) < 0.0 )
	{
	    this->_RotYY = acos(cosC2x);
        } else {
	    this->_RotYY = -acos(cosC2x);
	}
	LOG(BF("this->_RotYY = %lf degrees") % this->_RotYY/DEGTORAD  );
        this->_Distance = trans.length();
	LOG(BF("this->_Distance = %lf") % this->_Distance  );
#if 0 //[
//Vector3		translate;
//Vector3		labZ, labX, node, nodeNorm;
//double		cosRotY, gamma, beta, alpha, distance, theta, phi;
//	// Calculate RotZ, the angle around the Z-axis I need to rotate X into the lab XY-plane
//
//    labZ.set(0.0,0.0,1.0);
//    labX.set(1.0,0.0,0.0);
//    node = labZ.crossProduct(coord->getZ());
//    LOG(BF("labZ {%lf,%lf,%lf}") % labZ.getX() % labZ.getY() % labZ.getZ()  );
//    LOG(BF("labX {%lf,%lf,%lf}") % labX.getX() % labX.getY() % labX.getZ()  );
//    LOG(BF("coordX = {%lf,%lf,%lf}") % coord->getX().getX() % coord->getX().getY() % coord->getX().getZ()  );
//    LOG(BF("coordY = {%lf,%lf,%lf}") % coord->getY().getX() % coord->getY().getY() % coord->getY().getZ()  );
//    LOG(BF("coordZ = {%lf,%lf,%lf}") % coord->getZ().getX() % coord->getZ().getY() % coord->getZ().getZ()  );
//    cosRotY = labZ.dotProduct(coord->getZ());
//    LOG(BF("cosRotY = %lf") % cosRotY  );
//    LOG(BF("node.length= %lf") % node.length() );
//    if ( fabs(node.length())<PARALLEL_CUTOFF ) {
//	LOG(BF("node.length is zero") );
//	gamma = 0.0; // If the labZ and coord->getZ() are parallel or anti-parallel then gamma=
//        LOG(BF("gamma = %lf deg") % gamma/DEGTORAD  );
//	beta = (cosRotY>=0.0)?0.0:(180.0*DEGTORAD);
//        LOG(BF("beta = %lf deg") % (beta/DEGTORAD ) );
//        LOG(BF("Calculating alpha") );
//        alpha = labX.angleToVectorAboutNormal(coord->getX(),labZ);
//        LOG(BF("alpha= %lf deg") % alpha/DEGTORAD  );
//    } else {
//	LOG(BF("node.length is not zero") );
//		// RotY range is [0,PI], arc
//	nodeNorm = node.normalized(_lisp);
//        vp0(("nodeNorm =     {%lf,%lf,%lf};",nodeNorm.getX(), nodeNorm.getY(), nodeNorm.getZ() ));
//	LOG(BF("About to calculate gamma") );
//	gamma = nodeNorm.angleToVectorAboutNormal(coord->getX(),coord->getZ());
//        LOG(BF("gamma = %lf deg") % gamma/DEGTORAD  );
//	beta = acos(cosRotY);
//        LOG(BF("beta = %lf deg") % beta/DEGTORAD  );
//        LOG(BF("Calculating alpha") );
//        alpha = labX.angleToVectorAboutNormal(nodeNorm,labZ);
//        LOG(BF("alpha= %lf deg") % alpha/DEGTORAD  );
//    }
//
//	// Now calculate RotYY and RotZZ
//
//    translate = coord->getOrigin();
//    LOG(BF("translate {%lf,%lf,%lf}") % translate.getX() % translate.getY() % translate.getZ()  );
//    distance = translate.length();
//    LOG(BF("distance = %lf") % distance  );
//    if ( distance < 0.00000001 ) {
//        LOG(BF("distance is zero") );
//	phi = 0.0;
//	theta = 0.0;
//    } else {
//        LOG(BF("distance is not zero") );
//	Vector3 tn = translate.normalized(_lisp);
//	double cosRotZZ = tn.dotProduct(labZ);
//	LOG(BF("cosRotZZ = %lf") % cosRotZZ  );
//	Vector3 tnorm= labZ.crossProduct(tn);
//	if ( fabs(tnorm.length()) < PARALLEL_CUTOFF ) {
//	    LOG(BF("tnorm.length is zero ") );
//	    theta = (cosRotZZ>=0.0)?0.0:(180.0*DEGTORAD);
//	    phi = 0.0;
//	} else {
//	    LOG(BF("tnorm.length is not zero ") );
//	    theta = acos(cosRotZZ);
//	    phi = labX.angleToVectorAboutNormal(tnorm,labZ);
//	}
//	LOG(BF("theta = %lf deg") % theta/DEGTORAD  );
//	LOG(BF("phi = %lf deg") % phi/DEGTORAD  );
//    }
//    LOG(BF("Setting values") );
//    this->_RotX = alpha;
//    this->_RotY = beta;
//    this->_RotZ = gamma;
//    this->_Distance = distance;
//    this->_RotZZ = theta;
//    this->_RotYY = phi;
#endif //]
}

string	O_GimbalTransform::asString(const string& prefix)
{
stringstream	ss;
ss << prefix << "RotX  = " << this->_RotX/0.0174533 << " deg" << std::endl;
ss << prefix << "RotY  = " << this->_RotY/0.0174533 << " deg" << std::endl;
ss << prefix << "RotZ  = " << this->_RotZ/0.0174533 << " deg" << std::endl;
ss << prefix << "RotYY = " << this->_RotYY/0.0174533 << " deg" << std::endl;
ss << prefix << "RotZZ = " << this->_RotZZ/0.0174533 << " deg" << std::endl;
ss << prefix << "Dist  = " << this->_Distance << " A" << std::endl;
    return ss.str();
}

string	O_GimbalTransform::asLine(bool data)
{
stringstream	ss;
    if ( !data ) {
	ss << "# RotX(deg) RotY(deg) RotZ(deg) RotYY(deg) RotZZ(deg) Distance(A)";
    } else {
	ss << this->_RotX/0.0174533 << " ";
	ss << this->_RotY/0.0174533 << " " ;
	ss << this->_RotZ/0.0174533 << " " ;
	ss << this->_RotYY/0.0174533 << " ";
	ss << this->_RotZZ/0.0174533 << " ";
	ss << this->_Distance;
    }
    return ss.str();
}

void	O_GimbalTransform::defineForCoordinateSystemToCoordinateSystem(CoordinateSystem_sp start, CoordinateSystem_sp dest)
{
Matrix	toCanonical;
CoordinateSystem_sp	transformedDest;

    toCanonical = start->matrixToCanonical();
    transformedDest = dest->copy();
    transformedDest->transformWithMatrix(toCanonical);
    this->defineForCoordinateSystem(transformedDest);
}

void	O_GimbalTransform::defineFromGimbal(RPGimbalTransform g)
{
    this->_RotX = g->getRotX();
    this->_RotY = g->getRotY();
    this->_RotZ = g->getRotZ();
    this->_RotYY = g->getRotYY();
    this->_RotZZ = g->getRotZZ();
    this->_Distance = g->getDistance();
}

void	O_GimbalTransform::archive(core::ArchiveP n)
{
    n->attribute("RotX",this->_RotX );
    n->attribute("RotY",this->_RotY );
    n->attribute("RotZ",this->_RotZ );
    n->attribute("Dist",this->_Distance );
    n->attribute("RotYY",this->_RotYY );
    n->attribute("RotZZ",this->_RotZZ );
}

CoordinateSystem_sp	O_GimbalTransform::getCoordinateSystem()
{
Matrix		mt;
CoordinateSystem_sp	c;

c = CoordinateSystem_O::create();
    mt = this->matrixFromCanonical();
    c->transformWithMatrix(mt);
    return c;
}


    Vector3	O_GimbalTransform::getDirection()
{
    Vector3 vTrans; vTrans.set(1.0,0.0,0.0);
    Matrix mRotYY; mRotYY.rightHandedRotationY(this->_RotYY);
    Matrix mRotZZ; mRotZZ.rightHandedRotationZ(this->_RotZZ);
    vTrans = mRotYY.multiplyByVector3(vTrans);
    vTrans = mRotZZ.multiplyByVector3(vTrans);
    return vTrans;
}




Matrix	O_GimbalTransform::matrixFromCanonical()
{_G();
    Matrix mX; mX.rightHandedRotationX(this->_RotX);
    Matrix mY; mY.rightHandedRotationY(this->_RotY);
    Matrix mZ; mZ.rightHandedRotationZ(this->_RotZ);
    Vector3 vTrans = this->getDirection();
    vTrans = vTrans.multiplyByScalar(this->_Distance);
    Matrix mTranslate; mTranslate.translate(&vTrans);
    Matrix mA = mY*mX;
    mA = mZ*mA;
    mA = mTranslate*mA;
    return mA;
}


Matrix	O_GimbalTransform::matrixToCanonical()
{
CoordinateSystem_sp	c;
    c = this->getCoordinateSystem();
    return c->matrixToCanonical();
}





void	O_XConventionEulerTransform::defineForCoordinateSystem(CoordinateSystem_sp coord)
{_G();
    Vector3		translate;
    Vector3		labZ, labX, node, nodeNorm;
double		cosBeta, gamma, beta, alpha, distance, theta, phi;
	// Calculate Gamma, the angle around the Z-axis I need to rotate X into the lab XY-plane

    labZ.set(0.0,0.0,1.0);
    labX.set(1.0,0.0,0.0);
    node = labZ.crossProduct(coord->getZ());
    LOG(BF("labZ {%lf,%lf,%lf}") % labZ.getX() % labZ.getY() % labZ.getZ()  );
    LOG(BF("labX {%lf,%lf,%lf}") % labX.getX() % labX.getY() % labX.getZ()  );
    LOG(BF("coordX = {%lf,%lf,%lf}") % (coord->getX().getX()) % (coord->getX().getY()) % (coord->getX().getZ() ) );
    LOG(BF("coordY = {%lf,%lf,%lf}") % (coord->getY().getX()) % (coord->getY().getY()) % (coord->getY().getZ() ) );
    LOG(BF("coordZ = {%lf,%lf,%lf}") % (coord->getZ().getX()) % (coord->getZ().getY()) % (coord->getZ().getZ() ) );
    cosBeta = labZ.dotProduct(coord->getZ());
    LOG(BF("cosBeta = %lf") % cosBeta  );
    LOG(BF("node.length= %lf") % node.length() );
    if ( fabs(node.length())<PARALLEL_CUTOFF ) {
	LOG(BF("node.length is zero") );
	gamma = 0.0; // If the labZ and coord->getZ() are parallel or anti-parallel then gamma=
        LOG(BF("gamma = %lf deg") % gamma/DEGTORAD  );
	beta = (cosBeta>=0.0)?0.0:(180.0*DEGTORAD);
        LOG(BF("beta = %lf deg") % beta/DEGTORAD  );
        LOG(BF("Calculating alpha") );
        alpha = labX.angleToVectorAboutNormal(coord->getX(),labZ);
        LOG(BF("alpha= %lf deg") % alpha/DEGTORAD  );
    } else {
	LOG(BF("node.length is not zero") );
		// Beta range is [0,PI], arc
	nodeNorm = node.normalized(_lisp);
        LOG(BF("nodeNorm =     {%lf,%lf,%lf}") % (nodeNorm.getX()) % (nodeNorm.getY()) % (nodeNorm.getZ() ) );
	LOG(BF("About to calculate gamma") );
	gamma = nodeNorm.angleToVectorAboutNormal(coord->getX(),coord->getZ());
        LOG(BF("gamma = %lf deg") % gamma/DEGTORAD  );
	beta = acos(cosBeta);
        LOG(BF("beta = %lf deg") % beta/DEGTORAD  );
        LOG(BF("Calculating alpha") );
        alpha = labX.angleToVectorAboutNormal(nodeNorm,labZ);
        LOG(BF("alpha= %lf deg") % alpha/DEGTORAD  );
    }

	// Now calculate Phi and Theta

    translate = coord->getOrigin();
    LOG(BF("translate {%lf,%lf,%lf}") % translate.getX() % translate.getY() % translate.getZ()  );
    distance = translate.length();
    LOG(BF("distance = %lf") % distance  );
    if ( distance < 0.00000001 ) {
        LOG(BF("distance is zero") );
	phi = 0.0;
	theta = 0.0;
    } else {
        LOG(BF("distance is not zero") );
	Vector3 tn = translate.normalized(_lisp);
	double cosTheta = tn.dotProduct(labZ);
	LOG(BF("cosTheta = %lf") % cosTheta  );
	Vector3 tnorm= labZ.crossProduct(tn);
	if ( fabs(tnorm.length()) < PARALLEL_CUTOFF ) {
	    LOG(BF("tnorm.length is zero ") );
	    theta = (cosTheta>=0.0)?0.0:(180.0*DEGTORAD);
	    phi = 0.0;
	} else {
	    LOG(BF("tnorm.length is not zero ") );
	    theta = acos(cosTheta);
	    phi = labX.angleToVectorAboutNormal(tnorm,labZ);
	}
	LOG(BF("theta = %lf deg") % theta/DEGTORAD  );
	LOG(BF("phi = %lf deg") % phi/DEGTORAD  );
    }
    LOG(BF("Setting values") );
    this->_Alpha = alpha;
    this->_Beta = beta;
    this->_Gamma = gamma;
    this->_Distance = distance;
    this->_Theta = theta;
    this->_Phi = phi;
}

string	O_XConventionEulerTransform::asString()
{
stringstream	ss;
ss << "Alpha = " << this->_Alpha/0.0174533 << std::endl;
ss << "Beta = " << this->_Beta/0.0174533 << std::endl;
ss << "Gamma = " << this->_Gamma/0.0174533 << std::endl;
ss << "Phi = " << this->_Phi/0.0174533 << std::endl;
ss << "Theta = " << this->_Theta/0.0174533 << std::endl;
ss << "Distance = " << this->_Distance << std::endl;
    return ss.str();
}

void	O_XConventionEulerTransform::defineForCoordinateSystemToCoordinateSystem(CoordinateSystem_sp start, CoordinateSystem_sp dest)
{
Matrix	toCanonical;
CoordinateSystem_sp	transformedDest;

    toCanonical = start->matrixToCanonical();
    transformedDest = dest->copy();
    transformedDest->transformWithMatrix(toCanonical);
    this->defineForCoordinateSystem(transformedDest);
}

void	O_XConventionEulerTransform::archive(core::ArchiveP n)
{
    n->attribute("Alpha",this->_Alpha );
    n->attribute("Beta",this->_Beta );
    n->attribute("Gamma",this->_Gamma );
    n->attribute("Dist",this->_Distance );
    n->attribute("Phi",this->_Phi );
    n->attribute("Theta",this->_Theta );
}



CoordinateSystem_sp	O_XConventionEulerTransform::getCoordinateSystem()
{
Matrix		mt;
CoordinateSystem_sp	c;

c = CoordinateSystem_O::create();
    mt = this->matrixFromCanonical();
    c->transformWithMatrix(mt);
    return c;
}


Matrix	O_XConventionEulerTransform::matrixFromCanonical()
{_G();
    Matrix mD; mD.rightHandedRotationZ(this->_Alpha);
    Matrix mC; mC.rightHandedRotationX(this->_Beta);
    Matrix mB; mB.rightHandedRotationZ(this->_Gamma);
    Matrix mPhi; mPhi.rightHandedRotationZ(this->_Phi);
    Matrix mTheta; mTheta.rightHandedRotationX(this->_Theta);
    Vector3 vTrans; vTrans.set(0.0,0.0,this->_Distance);
    vTrans = mTheta.multiplyByVector3(vTrans);
    vTrans = mPhi.multiplyByVector3(vTrans);
    Matrix mTranslate; mTranslate.translate(&vTrans);
    Matrix mA = mC*mB;
    mA = mD*mA;
    mA = mTranslate*mA;
    return mA;
}


Matrix	O_XConventionEulerTransform::matrixToCanonical()
{
CoordinateSystem_sp	c;
    c = this->getCoordinateSystem();
    return c->matrixToCanonical();
}

#endif




void	CoordinateSystem_O::initialize()
{
    this->Base::initialize();
    this->canonical();
}


CoordinateSystem_O::CoordinateSystem_O(const CoordinateSystem_O& orig) : core::T_O(orig)
{
    this->origin = orig.origin;
    this->x = orig.x;
    this->y = orig.y;
    this->z = orig.z;
}


//
//	matrixFromCanonical
//
//	Generate the matrix for transformation
//	from the canonical coordinate system to this one
//
//	This means that a point at 0,0,0 will be transformed to
//	this coordinate systems origin.
//	A point at 1,0,0 will be transformed to a point
//	one unit along this ones X-axis
//
Matrix	CoordinateSystem_O::matrixFromCanonical()
{
Matrix	mt, m;

    mt.atColRowPut( 0, 0, this->x.getX() );
    mt.atColRowPut( 0, 1, this->x.getY() );
    mt.atColRowPut( 0, 2, this->x.getZ() );

    mt.atColRowPut( 1, 0, this->y.getX() );
    mt.atColRowPut( 1, 1, this->y.getY() );
    mt.atColRowPut( 1, 2, this->y.getZ() );

    mt.atColRowPut( 2, 0, this->z.getX() );
    mt.atColRowPut( 2, 1, this->z.getY() );
    mt.atColRowPut( 2, 2, this->z.getZ() );

    m.translate(&(this->origin));
    m = m*mt;
    return m;
}


//
//	matrixToCanonical
//
//	Generate the matrix for transformation
//	to the canonical coordinate system from this one
//
//	This means that a point on this coordinate systems origin will be
//	transformed to 0,0,0.  A point one unit along this coordinate systems
//	X-axis will be transformed to 1,0,0
//
Matrix	CoordinateSystem_O::matrixToCanonical()
{
Matrix  mt, m, mn;
Vector3 vt;
    m.atRowColPut( 0, 0, this->x.getX() );
    m.atRowColPut( 0, 1, this->x.getY() );
    m.atRowColPut( 0, 2, this->x.getZ() );

    m.atRowColPut( 1, 0, this->y.getX() );
    m.atRowColPut( 1, 1, this->y.getY() );
    m.atRowColPut( 1, 2, this->y.getZ() );

    m.atRowColPut( 2, 0, this->z.getX() );
    m.atRowColPut( 2, 1, this->z.getY() );
    m.atRowColPut( 2, 2, this->z.getZ() );

    vt = this->origin.multiplyByScalar(-1.0);

    mt.translate(&vt);
    LOG(BF("mt = %s") % (mt.asString().c_str() ) );
    mn = m*mt;
    LOG(BF("mn = %s") % (mn.asString().c_str() ) );
    return mn;
}


//
//	canonical
//
//	Set (this) to the canonical coordinate set.
//	origin is 0,0,0
//	x is 1,0,0
//	y is 0,1,0
//	z is 0,0,1
void CoordinateSystem_O::canonical()
{
    this->origin = Vector3(0.0,0.0,0.0);
    this->x = Vector3(1.0,0.0,0.0);
    this->y = Vector3(0.0,1.0,0.0);
    this->z = Vector3(0.0,0.0,1.0);
}


//
//	defineForAtoms
//
//	Build a coordinate set using the three atoms
//	the origin is on atom (aorigin)
//	the x axis is along (ax-aorigin)
//	the z axis is orthogonal to (axy) and x axis
//	the y axis is orthogonal to the new x and z axes
void	CoordinateSystem_O::defineForAtoms( Atom_sp aorigin, Atom_sp ax, Atom_sp axy )
{_G();
    Vector3	vo, vx, vxy, vz, vy;
    vo = aorigin->getPosition();
    LOG(BF("vo=%lf,%lf,%lf") % vo.getX() % vo.getY() % vo.getZ() );
    vx = (ax->getPosition() - vo).normalized();
    LOG(BF("vx=%lf,%lf,%lf") % vx.getX() % vx.getY() % vx.getZ() );
    vxy = (axy->getPosition() - vo).normalized();
    LOG(BF("vxy=%lf,%lf,%lf") % vxy.getX() % vxy.getY() % vxy.getZ() );
    vz = (vx.crossProduct(vxy)).normalized();
    LOG(BF("vz=%lf,%lf,%lf") % vz.getX() % vz.getY() % vz.getZ() );
    vy = (vz.crossProduct(vx)).normalized();
    LOG(BF("vy=%lf,%lf,%lf") % vy.getX() % vy.getY() % vy.getZ() );
    this->origin = vo;
    this->x = vx;
    this->y = vy;
    this->z = vz;
}


//
//	defineForAtomVectors
//
//	Build a coordinate set using the three atoms
//	the origin is on atom (aorigin)
//	the x axis is along (ax-aorigin)
//	the z axis is orthogonal to (axy) and x axis
//	the y axis is orthogonal to the new x and z axes
void	CoordinateSystem_O::defineForAtomVectors( const Vector3& vo,
						  const Vector3& vax,
						  const Vector3& vaxy )
{_G();
    Vector3	vx, vxy, vz, vy;
vx = (vax.sub(vo)).normalized();
vxy = (vaxy.sub(vo)).normalized();
    LOG(BF("vo=%lf,%lf,%lf") % vo.getX() % vo.getY() % vo.getZ() );
    LOG(BF("vx=%lf,%lf,%lf") % vx.getX() % vx.getY() % vx.getZ() );
    LOG(BF("vxy=%lf,%lf,%lf") % vxy.getX() % vxy.getY() % vxy.getZ() );

    vz = (vx.crossProduct(vxy)).normalized();
    LOG(BF("vz=%lf,%lf,%lf") % vz.getX() % vz.getY() % vz.getZ() );
    vy = (vz.crossProduct(vx)).normalized();
    LOG(BF("vy=%lf,%lf,%lf") % vy.getX() % vy.getY() % vy.getZ() );
    this->origin = vo;
    this->x = vx;
    this->y = vy;
    this->z = vz;
}

//
//	defineForTwoAtomVectors
//
//	Build a coordinate set using the two atoms
//	the origin is on atom (aorigin)
//	the x axis is along (ax-aorigin)
//	the y axis is whereever it ends up
//	the z axis is orthogonal to y and x axis
void	CoordinateSystem_O::defineForTwoAtomVectors( const Vector3& vo,
						     const Vector3& vax )
{_G();
    Vector3	vx, vz, vy;
vx = (vax.sub(vo)).normalized();
    LOG(BF("vo=%lf,%lf,%lf") % vo.getX() % vo.getY() % vo.getZ() );
    LOG(BF("vx=%lf,%lf,%lf") % vx.getX() % vx.getY() % vx.getZ() );
    Vector3 vOffset;
    vOffset.set(0.0,0.5,0.0);
    vy = (vx.add(vOffset)).normalized();
    	// Check if vy is a little offset from vx
    vz = ( vx.crossProduct(vy) );
    if ( vz.dotProduct(vz) < 0.0001 )
    {
	vOffset.set(0.5,0.0,0.0);
	vy = (vx.add(vOffset)).normalized();
	vz = ( vx.crossProduct(vy) ).normalized();
    } else
    {
	vz = vz.normalized();
    }
    vy = (vz.crossProduct(vx)).normalized();
    LOG(BF("vy=%lf,%lf,%lf") % vy.getX() % vy.getY() % vy.getZ() );
    LOG(BF("vz=%lf,%lf,%lf") % vz.getX() % vz.getY() % vz.getZ() );
    this->origin = vo;
    this->x = vx;
    this->y = vy;
    this->z = vz;
}

//
//	defineForOneAtomVectors
//
//	Build a coordinate set using the two atoms
//	the origin is on atom (aorigin)
//	the x axis is along (ax-aorigin)
//	the y axis is whereever it ends up
//	the z axis is orthogonal to y and x axis
void	CoordinateSystem_O::defineForOneAtomVector( const Vector3& vo )
{_G();
    Vector3	vx, vz, vy;
    this->origin = vo;
    this->x.set(1.0,0.0,0.0);
    this->y.set(0.0,1.0,0.0);
    this->z.set(0.0,0.0,1.0);
}



//
//	defineForVectorsOriginXDirYDir
//
//	Build a coordinate set using the origin and the direction
//	of the X-axis and Y-axis.  The Z-axis is calculated from
//	the cross-product of xDir and yDir.
//
void	CoordinateSystem_O::defineForVectorsOriginXDirYDir(
    const Vector3& orig,
    const Vector3& xDir,
    const Vector3& yDir )
{
    Vector3	vzt;
    this->origin = orig;
    this->x = xDir.normalized();
    this->y = yDir.normalized();
    vzt = (this->x).crossProduct(this->y);
    this->z = vzt.normalized();
}


//
//	defineForVectorsOriginXDirXYPlane
//
//	Build a coordinate set using the origin and the direction
//	of the X-axis and a vector in the XY plane.
//	The Z-axis is calculated from the cross-product of xDir and xyPlane.
//
void	CoordinateSystem_O::defineForVectorsOriginXDirXYPlane(
    const Vector3& orig,
    const Vector3& xDir,
    const Vector3& xyPlane )
{
    Vector3	vzt,vxy,vyt;
    this->origin = orig;
    this->x = xDir.normalized();
    vxy = xyPlane.normalized();
    vzt = (this->x).crossProduct(vxy);
    this->z = vzt.normalized();
    vyt = (this->z).crossProduct(this->x);
    this->y = vyt.normalized();
}

//
//	defineForVectorsOriginZDirXZPlane
//
//	Build a coordinate set using the origin and the direction
//	of the Z-axis and a vector in the XZ plane.
//	The Y-axis is calculated from the cross-product of xDir and xyPlane.
//
void	CoordinateSystem_O::defineForVectorsOriginZDirXZPlane(
    const Vector3& orig,
    const Vector3& zDir,
    const Vector3& xzPlane )
{
    Vector3		vyt,vxz,vxt;
    this->origin = orig;
    this->z = zDir.normalized();
    vxz = xzPlane.normalized();
    vyt = (this->z).crossProduct(vxz);
    this->y = vyt.normalized();
    vxt = (this->y).crossProduct(this->z);
    this->x = vxt.normalized();
}




//
//	defineForVectorsOriginXDirZDir
//
//	Build a coordinate set using the origin and the direction
//	of the X-axis and Z-axis.  The Y-axis is calculated from
//	the cross-product of xDir and zDir.
//
void	CoordinateSystem_O::defineForVectorsOriginXDirZDir(
    const Vector3& orig,
    const Vector3& xDir,
    const Vector3& zDir )
{
    Vector3	vyt;
    this->origin = orig;
    this->x = xDir.normalized();
    this->z = zDir.normalized();
    vyt = (this->z).crossProduct(this->x);
    this->y = vyt.normalized();
}


void	CoordinateSystem_O::defineRandom( )
{_OF();
double x,y,z;
#define	RANGE	20.0
x = core::randomNumber01()*RANGE - (RANGE/2.0);
y = core::randomNumber01()*RANGE - (RANGE/2.0);
z = core::randomNumber01()*RANGE - (RANGE/2.0);
    this->origin.set(x,y,z);
    while (1) {
	x = core::randomNumber01()*2.0-1.0;
	y = core::randomNumber01()*2.0-1.0;
	z = core::randomNumber01()*2.0-1.0;
	this->x.set(x,y,z);
	if ( this->x.length() > 0.00001 ) break;
    }
    this->x = this->x.normalized();
    while (1) {
	x = core::randomNumber01()*2.0-1.0;
	y = core::randomNumber01()*2.0-1.0;
	z = core::randomNumber01()*2.0-1.0;
	this->y.set(x,y,z);
	if ( this->y.length() > 0.00001 ) break;
    }
    this->z = this->x.crossProduct(this->y);
    this->z = this->z.normalized();
    this->y = this->z.crossProduct(this->x);
}


void	CoordinateSystem_O::defineFromCoordinateSystem( CoordinateSystem_sp c )
{
    this->origin = c->getOrigin();
    this->x = c->getX();
    this->y = c->getY();
    this->z = c->getZ();
}


void	CoordinateSystem_O::defineCanonical( )
{
    this->origin.set(0.0,0.0,0.0);
    this->x.set(1.0,0.0,0.0);
    this->y.set(0.0,1.0,0.0);
    this->z.set(0.0,0.0,1.0);
}








//
//	matrixForTransformTo
//
//	Return the matrix to transform from (this) to (cs)
Matrix	CoordinateSystem_O::matrixForTransformTo(CoordinateSystem_sp cs )
{
Matrix	m1, m2;

    m1 = this->matrixToCanonical();
    m2 = cs->matrixFromCanonical();
    return(m1*m2);
}




//
//	transformWithMatrix
//
//	Transform this coordinate system with the matrix
//
void	CoordinateSystem_O::transformWithMatrix(const Matrix& m )
{
    Vector3	vox, voy, voz;

    vox = m*(this->x + this->origin);
    voy = m*(this->y + this->origin);
    voz = m*(this->z + this->origin);
    this->origin = m*this->origin;
    this->x = vox - this->origin;
    this->y = voy - this->origin;
    this->z = voz - this->origin;
}

//
//	transformWithMatrixPointer
//
//	Transform this coordinate system with the matrix
//
void	CoordinateSystem_O::transformWithMatrixPointer( Matrix* m )
{
    Vector3	vox, voy, voz;

    vox = (*m)*(this->x + this->origin);
    voy = (*m)*(this->y + this->origin);
    voz = (*m)*(this->z + this->origin);
    this->origin = (*m)*this->origin;
    this->x = vox - this->origin;
    this->y = voy - this->origin;
    this->z = voz - this->origin;
}

#if 0
//
//	renderXml
//
//	Render the coordinate system in XML format
//
adapt::QDomNode_sp	CoordinateSystem_O::renderXml( double width,
							double axisLength)
{
IMPLEMENT_ME();
adapt::QDomNode_sp	graalphacs, line;
#if 0 //[

    graalphacs = xmlBlock();
    line = xmlLine( XML_RED, (int)width, this->origin,
				this->origin+this->x*axisLength );
    graalphacs->addChild(line);
    line = xmlLine( XML_GREEN, (int)width, this->origin,
				this->origin+this->y*axisLength );
    graalphacs->addChild(line);
    line = xmlLine( XML_BLUE, (int)width, this->origin,
				this->origin+this->z*axisLength );
    graalphacs->addChild(line);
#endif //]
    return graalphacs;
}
#endif



//
//	dump
//
//	Dump the coordinate system
//
void	CoordinateSystem_O::dump()
{
    _lisp->print(BF("origin: ( %lf, %lf, %lf )") %	this->origin.getX()% this->origin.getY() % this->origin.getZ() );
    _lisp->print(BF("x-axis: ( %lf, %lf, %lf )") %	this->x.getX()% this->x.getY() % this->x.getZ() );
    _lisp->print(BF("y-axis: ( %lf, %lf, %lf )") %	this->y.getX()% this->y.getY() % this->y.getZ() );
    _lisp->print(BF("z-axis: ( %lf, %lf, %lf )") %	this->z.getX()% this->z.getY() % this->z.getZ() );
}

#ifdef XML_ARCHIVE
void	CoordinateSystem_O::archive(core::ArchiveP node)
{
    node->archivePlainObject<Vector3>( "origin", "Vector", this->origin );
    node->archivePlainObject<Vector3>( "x", "Vector", this->x );
    node->archivePlainObject<Vector3>( "y", "Vector", this->y );
    node->archivePlainObject<Vector3>( "z", "Vector", this->z );
}
#endif


CoordinateSystem_sp	CoordinateSystem_O::copy()
{
    GC_COPY(CoordinateSystem_O, n , *this); // = RP_Copy<CoordinateSystem_O>(this);
    return n;
}


bool	CoordinateSystem_O::sameAs(CoordinateSystem_sp c)
{
    if ( !this->origin.sameAs(c->getOrigin()) ) return false;
    if ( !this->x.sameAs(c->getX()) ) return false;
    if ( !this->y.sameAs(c->getY()) ) return false;
    if ( !this->z.sameAs(c->getZ()) ) return false;
    return true;
}

#ifdef RENDER
geom::Render_sp CoordinateSystem_O::rendered(core::List_sp koptions)
{_G();
    geom::GrLines_sp lines = geom::GrLines_O::create();
    uint idxOrigin = lines->appendVertex(this->origin,geom::Color_O::white(_lisp));
    uint idxX = lines->appendVertex(this->origin+this->x,geom::Color_O::red(_lisp));
    uint idxY = lines->appendVertex(this->origin+this->y,geom::Color_O::green(_lisp));
    uint idxZ = lines->appendVertex(this->origin+this->z,geom::Color_O::blue(_lisp));
    lines->appendLine(idxOrigin,idxX);
    lines->appendLine(idxOrigin,idxY);
    lines->appendLine(idxOrigin,idxZ);
    geom::DisplayList_sp prims = geom::DisplayList_O::create();
    prims->append(lines);
    return prims;
}
#endif

string	CoordinateSystem_O::asString()
{
stringstream	ss;
ss << "Origin: " << this->origin.asString() << std::endl;
ss << "xdir: " << this->x.asString() << std::endl;
ss << "ydir: " << this->y.asString() << std::endl;
ss << "zdir: " << this->z.asString() << std::endl;
    return ss.str();
}

Vector3	CoordinateSystem_O::vectorRelativeToYou(const Vector3& pos)
{
    Vector3	newPos;
Matrix	transformToCanonical;
    transformToCanonical = this->matrixToCanonical();
    newPos = transformToCanonical.multiplyByVector3(pos);
    return newPos;
}



void CoordinateSystem_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<CoordinateSystem_O>()
	.def("getOrigin", &chem::CoordinateSystem_O::getOrigin)
	.def("setOrigin", &chem::CoordinateSystem_O::setOrigin)
	.def("getX", &chem::CoordinateSystem_O::getX)
	.def("getY", &chem::CoordinateSystem_O::getY)
	.def("getZ", &chem::CoordinateSystem_O::getZ)
	.def("canonical", &chem::CoordinateSystem_O::canonical)
	.def("defineForAtomVectors", &chem::CoordinateSystem_O::defineForAtomVectors)
	.def("defineForAtoms", &chem::CoordinateSystem_O::defineForAtoms)
	.def("defineForVectorsOriginXDirYDir", &chem::CoordinateSystem_O::defineForVectorsOriginXDirYDir)
	.def("defineForVectorsOriginXDirXYPlane", &chem::CoordinateSystem_O::defineForVectorsOriginXDirXYPlane)
	.def("defineForVectorsOriginZDirXZPlane", &chem::CoordinateSystem_O::defineForVectorsOriginZDirXZPlane)
	.def("defineForVectorsOriginXDirZDir", &chem::CoordinateSystem_O::defineForVectorsOriginXDirZDir)
	.def("defineFromCoordinateSystem", &chem::CoordinateSystem_O::defineFromCoordinateSystem)
	.def("defineCanonical", &chem::CoordinateSystem_O::defineCanonical)
	.def("defineRandom", &chem::CoordinateSystem_O::defineRandom)
	    .def("sameAs", &chem::CoordinateSystem_O::sameAs,"","","",false)
	.def("matrixFromCanonical", &chem::CoordinateSystem_O::matrixFromCanonical)
	.def("matrixToCanonical", &chem::CoordinateSystem_O::matrixToCanonical)
//	.def("renderXml", &chem::CoordinateSystem_O::renderXml)
	    .def("dump", &chem::CoordinateSystem_O::dump,"","","",false)
	    .def("copy", &chem::CoordinateSystem_O::copy,"","","",false)
//	.def("asXml", &chem::CoordinateSystem_O::asXml)
	.def("transformWithMatrix", &chem::CoordinateSystem_O::transformWithMatrix)
	;
    }

void CoordinateSystem_O::exposePython(core::Lisp_sp lisp)
{_G();
#ifdef USEBOOSTPYTHON
    PYTHON_CLASS(ChemPkg,CoordinateSystem,"","",_lisp)
	.def("getOrigin", &chem::CoordinateSystem_O::getOrigin)
	.def("setOrigin", &chem::CoordinateSystem_O::setOrigin)
	.def("getX", &chem::CoordinateSystem_O::getX)
	.def("getY", &chem::CoordinateSystem_O::getY)
	.def("getZ", &chem::CoordinateSystem_O::getZ)
	.def("canonical", &chem::CoordinateSystem_O::canonical)
	.def("defineForAtomVectors", &chem::CoordinateSystem_O::defineForAtomVectors)
	.def("defineForAtoms", &chem::CoordinateSystem_O::defineForAtoms)
	.def("defineForVectorsOriginXDirYDir", &chem::CoordinateSystem_O::defineForVectorsOriginXDirYDir)
	.def("defineForVectorsOriginXDirXYPlane", &chem::CoordinateSystem_O::defineForVectorsOriginXDirXYPlane)
	.def("defineForVectorsOriginZDirXZPlane", &chem::CoordinateSystem_O::defineForVectorsOriginZDirXZPlane)
	.def("defineForVectorsOriginXDirZDir", &chem::CoordinateSystem_O::defineForVectorsOriginXDirZDir)
	.def("defineFromCoordinateSystem", &chem::CoordinateSystem_O::defineFromCoordinateSystem)
	.def("defineCanonical", &chem::CoordinateSystem_O::defineCanonical)
	.def("defineRandom", &chem::CoordinateSystem_O::defineRandom)
//	.def("sameAs", &chem::CoordinateSystem_O::sameAs,"","","",false)
	.def("matrixFromCanonical", &chem::CoordinateSystem_O::matrixFromCanonical)
	.def("matrixToCanonical", &chem::CoordinateSystem_O::matrixToCanonical)
//	.def("renderXml", &chem::CoordinateSystem_O::renderXml)
//	.def("dump", &chem::CoordinateSystem_O::dump,"","","",false)
//	.def("copy", &chem::CoordinateSystem_O::copy,"","","",false)
//	.def("asXml", &chem::CoordinateSystem_O::asXml)
	.def("transformWithMatrix", &chem::CoordinateSystem_O::transformWithMatrix)
	;
#endif
    }


EXPOSE_CLASS(chem, CoordinateSystem_O);
}; // namespace chem


