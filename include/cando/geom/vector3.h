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

#include <limits>
#include <iostream>
#include <vector>
#include <clasp/core/foundation.h>
#include <clasp/core/lispStream.fwd.h>
#include <cando/chem/nVector.fwd.h>
//#include "core/archive.fwd.h"


#define MY_PI 3.1415926535897932
#define SMALL_NUMBER 1e-6

struct Safe {};
struct Unsafe {};

/*! A vector stores xyz coordinates in Angstroms */
class Vector3 {
  friend class OVector3_O;
protected:
  vecreal coords[3];
public:

  //! Construct zerod vector
  Vector3() : coords{
        std::numeric_limits<vecreal>::quiet_NaN(),
        std::numeric_limits<vecreal>::quiet_NaN(),
        std::numeric_limits<vecreal>::quiet_NaN() 
    } {};
  ~Vector3() {};
  Vector3( vecreal x, vecreal y, vecreal z ) : coords{x,y,z} {};

  Vector3( core::SimpleVector_float_sp svf, size_t index, Unsafe unsafe ) {
    ASSERT(index + 2 < svf->length());
    this->coords[0] = (*svf)[index];
    this->coords[1] = (*svf)[index + 1];
    this->coords[2] = (*svf)[index + 2];
  }
  Vector3( core::SimpleVector_float_sp svf, size_t index, Safe safe ) {
    if (index + 2 < svf->length()) {
      this->coords[0] = (*svf)[index];
      this->coords[1] = (*svf)[index + 1];
      this->coords[2] = (*svf)[index + 2];
      return;
    }
    SIMPLE_ERROR("The index {} is out of range for a 3D vector within the array of length {}", index, svf->length());
  }

  Vector3( core::SimpleVector_double_sp svf, size_t index, Unsafe unsafe ) {
    ASSERT(index + 2 < svf->length());
    this->coords[0] = (*svf)[index];
    this->coords[1] = (*svf)[index + 1];
    this->coords[2] = (*svf)[index + 2];
  }

  Vector3( core::SimpleVector_double_sp svf, size_t index, Safe safe ) {
    if (index + 2 < svf->length()) {
      this->coords[0] = (*svf)[index];
      this->coords[1] = (*svf)[index + 1];
      this->coords[2] = (*svf)[index + 2];
      return;
    }
    SIMPLE_ERROR("The index {} is out of range for a 3D vector within the array of length {}", index, svf->length());
  }


  /*! Return the vector in nanometers */
  Vector3 inNanometers() const;

  vecreal& operator[](uint x) { return this->coords[x];};

  vecreal const& operator[](uint x) const { return this->coords[x];};

  vecreal& getX() { return this->coords[0]; }
  vecreal& getY() { return this->coords[1]; }
  vecreal& getZ() { return this->coords[2]; }
  vecreal getX() const { return this->coords[0]; }
  vecreal getY() const { return this->coords[1]; }
  vecreal getZ() const { return this->coords[2]; }
  bool operator==( const Vector3& v ) const {
    if (std::isnan(this->coords[0])) {
      return std::isnan(v.coords[0]);
    } else {}
    return this->coords[0] == v.coords[0]
        && this->coords[1] == v.coords[1]
        && this->coords[2] == v.coords[2];
  }
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
  void	setIsDefined(bool isDef) { if ( !isDef ) { this->set(
          std::numeric_limits<vecreal>::quiet_NaN(),
          std::numeric_limits<vecreal>::quiet_NaN(),
          std::numeric_limits<vecreal>::quiet_NaN()
                                                             ); }; };
  //! A Vector3 is defined if any of its coordinates are not zero
  bool isDefined() const { return ( !std::isnan(this->coords[0]) ); }

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
  double	calculateDistance( const Vector3& va,
                                   const Vector3& vb );
  double	calculateDistanceSquared( const Vector3& va,
                                          const Vector3& vb);

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
