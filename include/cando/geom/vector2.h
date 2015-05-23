       
       
//
// (C) 2004 Christian E. Schafmeister
//



#ifndef	Vector2_H
#define Vector2_H

#include <iostream>
#include <vector>
#include <clasp/core/foundation.h>
//#include "core/archiveNode.h"


namespace geom {





class Vector2 {

public:

	Vector2();
	virtual ~Vector2();
	Vector2(double x, double y );
	double& operator[](uint x) { return this->coords[x];};

	Vector2 operator+( Vector2 v ) {
		return Vector2( this->coords[0]+v.getX(),
				 this->coords[1]+v.getY() );
	};
	Vector2 operator- ( Vector2 v ) const {
		return Vector2( this->coords[0]-v.getX(),
				 this->coords[1]-v.getY() );
	};
	Vector2 operator*( double s ) {
		return Vector2( this->coords[0]*s,
				this->coords[1]*s );
	};

	bool isDefined();
	void setIsDefined(bool b) {if (!b) {this->set(0.0,0.0);};};
	Vector2	add( Vector2 v ) { return (*this)+v; };
	Vector2	sub( Vector2 v ) { return (*this)-v; };
	Vector2	multiplyByScalar( double v ) { return (*this)*v; };

	bool isVector2() { return true; };
	void	set( double x, double y ) {
			this->coords[0]=x;
			this->coords[1]=y;
	};
	double	dotProduct( Vector2 v );
	double	length() const;
    Vector2	normalized(core::Lisp_sp lisp) const;

	void	dump();
//    void	archive( core::ArchiveP node );

	int	len() { return 2;};
	double get(int i) const;
	double getX() const { return this->coords[0]; }
	double getY() const { return this->coords[1]; }
	string	asString();


private:
	double coords[2];
};


typedef	vector<Vector2>	VectorVector2s;

double	calculateAngle( const Vector2& va,
			const Vector2& vb,
			const Vector2& vc,
			const core::Lisp_sp& lisp );

double	calculateDistance( const Vector2& va,
			   const Vector2& vb,
			   const core::Lisp_sp& lisp );



};


#endif
