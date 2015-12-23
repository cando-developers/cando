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
	static void lisp_initGlobals(core::Lisp_sp lisp);
    public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
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
	static Color_sp getOrDefineSystemColor(core::ArchiveP node);
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


};
TRANSLATE(geom::Color_O);
#endif //]
