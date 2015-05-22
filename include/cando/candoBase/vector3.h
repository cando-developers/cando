//
// (C) 2004 Christian E. Schafmeister
//



#ifndef	Vector3_H
#define Vector3_H

#include <iostream>
#include <vector>
#include <clasp/core/foundation.h>
#include <clasp/core/lispStream.fwd.h>
//#include "core/archive.fwd.h"
//#include "core/serialize.h"

/*! A vector stores xyz coordinates in Angstroms */
class Vector3 {
    friend class OVector3_O;

protected:
	double coords[3];
public:

	Vector3();
	virtual ~Vector3();
	Vector3(double x, double y, double z );

	/*! Return the vector in nanometers */
	Vector3 inNanometers() const;

	double& operator[](uint x) { return this->coords[x];};

	double const& operator[](uint x) const { return this->coords[x];};

	double getX() const { return this->coords[0]; }
	double getY() const { return this->coords[1]; }
	double getZ() const { return this->coords[2]; }
	Vector3 operator+( const Vector3& v ) const { return Vector3( this->coords[0]+v.getX(),
							 this->coords[1]+v.getY(),
							 this->coords[2]+v.getZ() ); };
	Vector3 operator- ( const Vector3& v ) const { return Vector3( this->coords[0]-v.getX(),
							 this->coords[1]-v.getY(),
							 this->coords[2]-v.getZ() ); };
	Vector3 operator*( double s ) { return Vector3( this->coords[0]*s,
							this->coords[1]*s,
							this->coords[2]*s ); };

	Vector3	add( const Vector3& v ) const { return (*this)+v; }  ;
	void addToX(const double& offset) { this->coords[0] += offset;};
	void addToY(const double& offset) { this->coords[1] += offset;};
	void addToZ(const double& offset) { this->coords[2] += offset;};
	Vector3	sub( const Vector3& v )  const{ return (*this)-v; } ;
	Vector3	multiplyByScalar( double v ) { return (*this)*v; };

	void	set( double x, double y, double z )	{ this->coords[0]=x;
							this->coords[1]=y;
							this->coords[2]=z;};
	void	setVector( double x, double y, double z) {this->set(x,y,z);};
	Vector3 crossProduct( const Vector3& v ) const;
	double	dotProduct( const Vector3& v );
	double	length() const;

		//! If the Vector3 is not defined then set all coordinates to 0.0
	void	setIsDefined(bool isDef) { if ( !isDef ) { this->set(0.0,0.0,0.0); }; };
		//! A Vector3 is defined if any of its coordinates are not zero
    bool	isDefined() { return ( this->coords[0]!=0.0 || this->coords[1]!=0.0 || this->coords[2]!=0.0 ); };
    void	archive( core::ArchiveP node );

	/*! Normalize the Vector3 and throw a hard exception if it cannot be */
	Vector3	normalized() const;

	Vector3	normalized(const core::Lisp_sp& lisp) const;
		//! Normalize the vector and if it is zero then return a zero vector
	Vector3	normalizedOrZero() const;

	double	angleToVectorAboutNormal(const Vector3& toVector, const Vector3& aboutNormal, const core::Lisp_sp& lisp);
#if 0 //[
	QDomNode_sp	asXml(string nm="V3");
	void		parseFromXml(QDomNode_sp node);
#endif //]
	void	dump(const core::Lisp_sp& lisp);
	bool	sameAs(const Vector3& x);

    void	writeToStream( std::ostream& out );
    void write(core::T_sp out);
	string	asString() const;
	void	parseFromString(string& str);
	void	parseFromCharacterArray(char* line);

	/*! Fill from a cons of three Numbers */
	void	fillFromCons(core::Cons_sp vals,core::Lisp_sp lisp);
};


typedef	vector<Vector3>	VectorVector3s;



double	calculateDihedral( const Vector3& va,
				const Vector3& vb,
				const Vector3& vc,
			   const Vector3& vd,
			   const core::Lisp_sp& lisp);

double	calculateAngle( const Vector3& va,
				const Vector3& vb,
			const Vector3& vc,
			const core::Lisp_sp& lisp );

double	calculateDistance( const Vector3& va,
			   const Vector3& vb,
			   const core::Lisp_sp& lisp );

/*! Faster than calculateDistance because it doesnt
 * require a square root
 */
double	calculateDistanceSquared( const Vector3& va,
				  const Vector3& vb);



Vector3	buildOrigin(const core::Lisp_sp& lisp);

Vector3	buildUsingBond( double distance, const Vector3& vb,
			const core::Lisp_sp& lisp);

//! Build a vector at distance from vb and angle from v
Vector3 buildUsingBondAngle( double distance, const Vector3& vb,
			     double angleRad, const Vector3& va,
			     const core::Lisp_sp& lisp);

Vector3 buildUsingBondAngleDihedral( double distance, const Vector3& vc,
					double angleRad, const Vector3& vb,
				     double dihedralRad, const Vector3& va );



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
	ASSERTF(pod!=NULL,BF("Did not find POD:Vector3 - found[%s]") % node->className());
	pod->resetCursorToHeadChild();
	core::Symbol_sp head;
	pod->loadObject(head);
	ASSERTF(head->isKeywordSymbol()&&head->symbolName()=="V3",BF("Expected Vector3->V3 as head but got[%s]") % head->symbolName() );
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
	core::Fixnum_sp iv = core::Fixnum_O::create((int)vals.size());
	pod->saveObject(iv);
	for (uint i=0; i<vals.size(); i++ ) pod->savePOD(vals[i]);
	return pod;
    }


    template <>
    inline void ContainerNode::loadPOD(vector<Vector3>& vals)
    {_OF();
	Node_p node = this->nextChild();
	PODNode_p pod = dynamic_cast<PODNode_p>(node);
	ASSERTF(pod!=NULL,BF("Did not find POD:vector<Vector3> - found[%s]") % node->className());
	pod->resetCursorToHeadChild();
	core::Symbol_sp head;
	pod->loadObject(head);
	ASSERTF(head->isKeywordSymbol()&&head->symbolName()=="VectorVector3",BF("Expected VectorVector3 as head but got[%s]") % head->symbolName() );
	int sz;
	pod->loadPOD(sz);
	vals.resize(sz);
	for (uint i=0; i<vals.size(); i++ ) pod->loadPOD(vals[i]);
    }
};
#endif // defined(OLD_SERIALIZE)


#endif
