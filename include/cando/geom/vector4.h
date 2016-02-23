/*
    File: vector4.h
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
	Vector4	normalized();

	void	dump();
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
