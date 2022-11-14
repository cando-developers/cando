/*
    File: vector3.h
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
//
// (C) 2004 Christian E. Schafmeister
//



#ifndef	Vector3_H
#define Vector3_H

#include <iostream>
#include <vector>
#include <clasp/core/foundation.h>
#include <clasp/core/lispStream.fwd.h>
#include <cando/chem/nVector.fwd.h>
//#include "core/archive.fwd.h"


#define MY_PI 3.1415926535897932
#define SMALL_NUMBER 1e-6

/*! A vector stores xyz coordinates in Angstroms */
class Vector3 {
  friend class OVector3_O;

protected:
  double coords[3];
public:


  //! Construct zerod vector
  Vector3() : coords{0.0,0.0,0.0} {};
  ~Vector3() {};
  Vector3( double x, double y, double z ) : coords{x,y,z} {};

  /*! Return the vector in nanometers */
  Vector3 inNanometers() const;

  double& operator[](uint x) { return this->coords[x];};

  double const& operator[](uint x) const { return this->coords[x];};

  double& getX() { return this->coords[0]; }
  double& getY() { return this->coords[1]; }
  double& getZ() { return this->coords[2]; }
  double getX() const { return this->coords[0]; }
  double getY() const { return this->coords[1]; }
  double getZ() const { return this->coords[2]; }
  bool operator==( const Vector3& v ) const { return this->coords[0] == v.coords[0]
        && this->coords[1] == v.coords[1]
        && this->coords[2] == v.coords[2]; }
  Vector3 operator+( const Vector3& v ) const { return Vector3( this->coords[0]+v.getX(),
                                                                this->coords[1]+v.getY(),
                                                                this->coords[2]+v.getZ() ); };
  Vector3 operator- ( const Vector3& v ) const { return Vector3( this->coords[0]-v.getX(),
                                                                 this->coords[1]-v.getY(),
                                                                 this->coords[2]-v.getZ() ); };
  Vector3 operator*( double s ) { return Vector3( this->coords[0]*s,
                                                  this->coords[1]*s,
                                                  this->coords[2]*s ); };

  Vector3 add( const Vector3& v ) const { return (*this)+v; }  ;
  void addSet( Vector3& result, const Vector3& v ) const { result = (*this)+v; }  ;
  void addToX(const double& offset) { this->coords[0] += offset;};
  void addToY(const double& offset) { this->coords[1] += offset;};
  void addToZ(const double& offset) { this->coords[2] += offset;};
  Vector3 sub( const Vector3& v ) const{ return (*this)-v; } ;
  void subSet( Vector3& result, const Vector3& v ) const{ result = (*this)-v; } ;
  Vector3	multiplyByScalar( double v ) { return (*this)*v; };
  void multiplyByScalarSet( Vector3& result, double s) {
    result.coords[0] = this->coords[0]*s;
    result.coords[1] = this->coords[1]*s;
    result.coords[2] = this->coords[2]*s;
  }
  

  void	set( double x, double y, double z )
  {
    this->coords[0]=x;
    this->coords[1]=y;
    this->coords[2]=z;
  };
  void	setVector( double x, double y, double z) {this->set(x,y,z);};
  Vector3 crossProduct( const Vector3& v ) const;
  void crossProductSet(const Vector3& v, Vector3& result ) const;
  double dotProduct( const Vector3& v ) const {
    return(this->getX()*v.getX()
           +this->getY()*v.getY()
           +this->getZ()*v.getZ());
  }
  double	length() const;
  double	lengthSquared() const {
    return this->dotProduct(*this);
  }

  //! If the Vector3 is not defined then set all coordinates to 0.0
  void	setIsDefined(bool isDef) { if ( !isDef ) { this->set(0.0,0.0,0.0); }; };
  //! A Vector3 is defined if any of its coordinates are not zero
  bool	isDefined() { return ( this->coords[0]!=0.0 || this->coords[1]!=0.0 || this->coords[2]!=0.0 ); };

  Vector3 normalized() const;
  void normalizedSet(Vector3& result);
  //! Normalize the vector and if it is zero then return a zero vector
  Vector3 normalizedOrZero() const;
  void normalizedOrZeroSet(Vector3& result) const;

  double	angleToVectorAboutNormal(const Vector3& toVector, const Vector3& aboutNormal);
#if 0 //[
  QDomNode_sp	asXml(string nm="V3");
  void		parseFromXml(QDomNode_sp node);
#endif //]
  void	dump();
  bool	sameAs(const Vector3& x);

  void	writeToStream( std::ostream& out );
  void write(core::T_sp out);
  string	asString() const;
  void	parseFromString(string& str);
  void	parseFromCharacterArray(char* line);

  /*! Fill from a cons of three Numbers */
  void	fillFromCons(core::Cons_sp vals);
};


typedef	vector<Vector3>	VectorVector3s;



/*! Faster than calculateDistance because it doesnt
 * require a square root
 */
namespace geom {
  double	calculateDihedral( const Vector3& va,
                                   const Vector3& vb,
                                   const Vector3& vc,
                                   const Vector3& vd );
  double	calculateAngle( const Vector3& va,
				const Vector3& vb,
                                const Vector3& vc );

  double	calculateDistance( const Vector3& va,
                                   const Vector3& vb );
  double	calculateDistanceSquared( const Vector3& va,
                                          const Vector3& vb);
Vector3	geom__build_origin();
  Vector3	geom__build_using_bond( double distance, const Vector3& vb );
//! Build a vector at distance from vb and angle from v
  Vector3 geom__build_using_bond_angle( double distance, const Vector3& vb,
                               double angleRad, const Vector3& va );
  Vector3 geom__build_using_bond_angle_dihedral( double distance, const Vector3& vc,
                                       double angleRad, const Vector3& vb,
                                       double dihedralRad, const Vector3& va );

/*! Return the angle in radians between (1.0,0.0) and (dx,dy). The result will be positive or negative depending on the quadrant of dx,dy */
double  geom__planeVectorAngle(double dx, double dy);

};

#if defined(OLD_SERIALIZE)
namespace serialize
{
// load/save vector<Vector3>'s
    template <>
    inline Node_p ContainerNode::savePOD(Vector3& vec)
    {_OF();
	PODNode_p pod = this->insertPODNodeAfterCursorAndAdvanceCursor();
	core::Symbol_sp head = _lisp->internKeyword("V3");
	pod->saveObject(head);
	CR_HINT(pod,false);
	pod->savePOD(vec[0]);
	CR_HINT(pod,false);
	pod->savePOD(vec[1]);
	CR_HINT(pod,false);
	pod->savePOD(vec[2]);
	CR_HINT(pod,false);
	return pod;
    }


    template <>
    inline void ContainerNode::loadPOD(Vector3& vec)
    {_OF();
	Node_p node = this->nextChild();
	PODNode_p pod = dynamic_cast<PODNode_p>(node);
	ASSERTF(pod!=NULL,("Did not find POD:Vector3 - found[%s]") , node->className());
	pod->resetCursorToHeadChild();
	core::Symbol_sp head;
	pod->loadObject(head);
	ASSERTF(head->isKeywordSymbol()&&head->symbolName()=="V3",("Expected Vector3->V3 as head but got[%s]") , head->symbolName() );
	pod->loadPOD(vec[0]);
	pod->loadPOD(vec[1]);
	pod->loadPOD(vec[2]);
    }




// load/save vector<Vector3>'s
    template <>
    inline Node_p ContainerNode::savePOD(vector<Vector3>& vals)
    {_OF();
	PODNode_p pod = this->insertPODNodeAfterCursorAndAdvanceCursor();
	core::Symbol_sp head = _lisp->internKeyword("VectorVector3");
	pod->saveObject(head);
	core::Fixnum_sp iv = core::clasp_make_fixnum(vals.size());
	pod->saveObject(iv);
	for (uint i=0; i<vals.size(); i++ ) pod->savePOD(vals[i]);
	return pod;
    }


    template <>
    inline void ContainerNode::loadPOD(vector<Vector3>& vals)
    {_OF();
	Node_p node = this->nextChild();
	PODNode_p pod = dynamic_cast<PODNode_p>(node);
	ASSERTF(pod!=NULL,("Did not find POD:vector<Vector3> - found[%s]") , node->className());
	pod->resetCursorToHeadChild();
	core::Symbol_sp head;
	pod->loadObject(head);
	ASSERTF(head->isKeywordSymbol()&&head->symbolName()=="VectorVector3",("Expected VectorVector3 as head but got[%s]") , head->symbolName() );
	int sz;
	pod->loadPOD(sz);
	vals.resize(sz);
	for (uint i=0; i<vals.size(); i++ ) pod->loadPOD(vals[i]);
    }
};
#endif // defined(OLD_SERIALIZE)

namespace geom {
void geom__vec_extract(Vector3& vec, chem::NVector_sp coordinates, size_t index0);

};

#endif
