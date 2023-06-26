/*
    File: vector2.h
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



#ifndef	Vector2_H
#define Vector2_H

#include <iostream>
#include <vector>
#include <clasp/core/foundation.h>
#include <cando/geom/vector3.fwd.h>
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
    Vector2	normalized() const;

	void	dump();
	int	len() { return 2;};
	double get(int i) const;
	double getX() const { return this->coords[0]; }
	double getY() const { return this->coords[1]; }
  string	asString() const;


private:
	double coords[2];
};


typedef	vector<Vector2>	VectorVector2s;
double	calculate2Angle( const Vector2& va,
			const Vector2& vb,
			const Vector2& vc);

double	calculate2Distance( const Vector2& va,
			   const Vector2& vb);


};


#endif
