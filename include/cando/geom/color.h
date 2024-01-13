/*
    File: color.h
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
#ifndef	Color_H //[
#define Color_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/color.fwd.h>

namespace geom
{


    SMART(Color);
    class Color_O : public core::CxxObject_O
    {
	LISP_CLASS(geom,GeomPkg,Color_O,"Color",core::CxxObject_O);
    public:
    public: // virtual functions inherited from Object
	void	initialize();
//	string	__repr__() const;

    private: // instance variables
	core::Symbol_sp	_ColorName;
	uint		_Color;
    public:
	/*! Construct with an RGB value */
	static Color_sp create(core::Symbol_sp name, uint rgb);

    public:
	static Color_sp systemColor(core::Symbol_sp name);
	static Color_sp getOrDefineSystemColor(core::Symbol_sp name, uint rgb);
	static Color_sp getOrDefineSystemColor(core::Symbol_sp name, double red, double green, double blue);
    public: // standard colors
	static Color_sp white();
	static Color_sp yellow();
	static Color_sp red();
	static Color_sp green();
	static Color_sp blue();
	static Color_sp cyan();
	static Color_sp magenta();
	static Color_sp black();
    public:
	void setColorName(core::Symbol_sp co);
	Color_sp setColor(uint c) { this->_Color = c; return this->sharedThis<Color_O>();};
	uint getColor() { return this->_Color;};
	string asString();
	void	setName(core::Symbol_sp sym);
	void	fillDoubleArray(double values[4]);
	void	fillFloatArray(float values[4]);
	void setRGB(double r, double g, double b);
	double	getRed();
	double	getGreen();
	double	getBlue();
	double	getAlpha();
	core::ubyte	getRedByte();
	core::ubyte	getGreenByte();
	core::ubyte	getBlueByte();
	core::ubyte	getAlphaByte();
	Color_O(const Color_O& orig);

    public:

	DEFAULT_CTOR_DTOR(Color_O);
    };


    void initialize_colors();

};
#endif //]
