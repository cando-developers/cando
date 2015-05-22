       
       
//
// (C) 2004 Christian E. Schafmeister
//



#ifndef	Vector4_H
#define Vector4_H

#include <iostream>
#include <vector>
#include <clasp/core/foundation.h>


class Vector4 {

public:

	Vector4();
	virtual ~Vector4();
	Vector4(double w, double x, double y, double z );
	Vector4 operator+( Vector4 v ) { return Vector4(
						this->coords[0]+v.getW(),
						this->coords[1]+v.getX(),
						this->coords[2]+v.getY(),
						this->coords[3]+v.getZ() ); };
	Vector4 operator-( Vector4 v ) { return Vector4(
						this->coords[0]-v.getW(),
						this->coords[1]-v.getX(),
						this->coords[2]-v.getY(),
						this->coords[3]-v.getZ() ); };
	Vector4 operator*( double s ) { return Vector4(
						this->coords[0]*s,
						this->coords[1]*s,
						this->coords[2]*s,
						this->coords[3]*s ); };

	Vector4	add( Vector4 v ) { return (*this)+v; };
	Vector4	sub( Vector4 v ) { return (*this)-v; };
	Vector4	multiplyByScalar( double v ) { return (*this)*v; };

	void	set( double w, double x, double y,double z )	{
							this->coords[0]=w;
							this->coords[1]=x;
							this->coords[2]=y;
							this->coords[3]=z;
	};
	double	dotProduct( Vector4 v );
	double	length();
	Vector4	normalized(const core::Lisp_sp& lisp);

	void	dump(const core::Lisp_sp& lisp);
	string	asString();

	double getW() { return this->coords[0]; }
	double getX() { return this->coords[1]; }
	double getY() { return this->coords[2]; }
	double getZ() { return this->coords[3]; }

	int	indexOfLargestElement();

    void	writeToStream( std::ostream& out );
private:
	double coords[4];
};


typedef	vector<Vector4>	VectorVector4s;



#endif
