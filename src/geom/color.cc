/*
    File: color.cc
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
#define	DEBUG_LEVEL_NONE

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/symbolTable.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/lisp.h>
#include <cando/geom/color.h>

// last include is wrappers.h
#include <clasp/core/wrappers.h>


namespace geom
{
    SYMBOL_EXPORT_SC_(GeomPkg,colorTable);
    SYMBOL_EXPORT_SC_(KeywordPkg,white);
    SYMBOL_EXPORT_SC_(KeywordPkg,black);
    SYMBOL_EXPORT_SC_(KeywordPkg,red);
    SYMBOL_EXPORT_SC_(KeywordPkg,green);
    SYMBOL_EXPORT_SC_(KeywordPkg,blue);
    SYMBOL_EXPORT_SC_(KeywordPkg,magenta);
    SYMBOL_EXPORT_SC_(KeywordPkg,yellow);
    SYMBOL_EXPORT_SC_(KeywordPkg,cyan);
    SYMBOL_EXPORT_SC_(KeywordPkg,grey10);
    SYMBOL_EXPORT_SC_(KeywordPkg,greyF0);
    SYMBOL_EXPORT_SC_(KeywordPkg,goldenrod);
    SYMBOL_EXPORT_SC_(KeywordPkg,purple);
    SYMBOL_EXPORT_SC_(KeywordPkg,pink);
    SYMBOL_EXPORT_SC_(KeywordPkg,plum);
    SYMBOL_EXPORT_SC_(KeywordPkg,coral);
    SYMBOL_EXPORT_SC_(KeywordPkg,orange);
    SYMBOL_EXPORT_SC_(KeywordPkg,brown);
    SYMBOL_EXPORT_SC_(KeywordPkg,wheat);
    SYMBOL_EXPORT_SC_(KeywordPkg,chartreuse);
    SYMBOL_EXPORT_SC_(KeywordPkg,aquamarine);
    SYMBOL_EXPORT_SC_(KeywordPkg,mediumBlue);
    SYMBOL_EXPORT_SC_(KeywordPkg,steelBlue);





void initialize_colors()
{
  core::HashTable_sp colorTable = core::HashTable_O::createEq();
  SYMBOL_EXPORT_SC_(GeomPkg,colorTable);
  _lisp->defvar(_sym_colorTable,colorTable);
#define DEF_COLOR(name,rgb) {\
    Color_O::getOrDefineSystemColor(kw::_sym_ ## name, rgb); \
  }
  DEF_COLOR(white,0xffFFff);
  DEF_COLOR(black,0x000000);
  DEF_COLOR(red,0xFF0000);
  DEF_COLOR(green,0x00FF00);
  DEF_COLOR(blue,0x0000FF);
  DEF_COLOR(magenta,0xFF00FF);
  DEF_COLOR(yellow,0xFFFF00);
  DEF_COLOR(cyan,0x00FFFF);
  DEF_COLOR(grey10,0x101010);
  DEF_COLOR(greyF0,0xF0F0F0);
  DEF_COLOR(goldenrod,0xdaa520);
  DEF_COLOR(purple,0xa020f0);
  DEF_COLOR(pink,0xffc0cb);
  DEF_COLOR(plum,0xdda0dd);
  DEF_COLOR(coral,0xff7f50);
  DEF_COLOR(orange,0xffa500);
  DEF_COLOR(brown,0xa52a2a);
  DEF_COLOR(wheat,0xf5deb3);
  DEF_COLOR(chartreuse,0x7ff000);
  DEF_COLOR(aquamarine,0x7fffd4);
  DEF_COLOR(mediumBlue,0x0000cd);
  DEF_COLOR(steelBlue,0x4682b4);
};





Color_O::Color_O(const Color_O& orig) : core::CxxObject_O(orig)
{
    this->_ColorName = orig._ColorName;
    this->_Color = orig._Color;
}



    Color_sp Color_O::white(){return Color_O::systemColor(kw::_sym_white);}
    Color_sp Color_O::yellow(){return Color_O::systemColor(kw::_sym_yellow);}
    Color_sp Color_O::red(){return Color_O::systemColor(kw::_sym_red);}
    Color_sp Color_O::green(){return Color_O::systemColor(kw::_sym_green);}
    Color_sp Color_O::blue(){return Color_O::systemColor(kw::_sym_blue);}
    Color_sp Color_O::cyan(){return Color_O::systemColor(kw::_sym_cyan);}
    Color_sp Color_O::magenta(){return Color_O::systemColor(kw::_sym_magenta);}
    Color_sp Color_O::black(){return Color_O::systemColor(kw::_sym_black);}


    Color_sp Color_O::systemColor(core::Symbol_sp sym)
    {
	core::Symbol_sp colorTableSymbol = _sym_colorTable;
	core::HashTable_sp colorTable = colorTableSymbol->symbolValue().as<core::HashTable_O>();
        core::T_sp c = colorTable->gethash(sym);
        if (c.nilp()) {
          SIMPLE_ERROR("Could not find color {}" , _rep_(sym));
        }
	return c.as<Color_O>();
    }


    Color_sp Color_O::getOrDefineSystemColor(core::Symbol_sp sym, uint color)
{
    core::HashTable_sp colorTable = _sym_colorTable->symbolValue().as<core::HashTable_O>();
    core::T_sp c = colorTable->gethash(sym);
    if ( c.nilp() ) {
	Color_sp newColor = Color_O::create(sym,color);
	colorTable->setf_gethash(sym,newColor);
	return newColor;
    }
    return c.as<Color_O>();
}
    
#ifdef XML_ARCHIVE
    Color_sp Color_O::getOrDefineSystemColor(core::ArchiveP node)
{
    core::HashTable_sp colorTable = _sym_colorTable->symbolValue().as<core::HashTable_O>();
    string colorNameString = node->getAttributeString("colorName");
    core::Symbol_sp sym = _lisp->internKeyword(colorNameString);
    core::Binder_O::iterator it = colorTable->find(sym);
    if ( it == colorTable->end() )
    {
	uint uicolor = node->getAttributeUInt("color");
	Color_sp newColor = Color_O::create(sym,uicolor);
	colorTable->extend(sym,newColor);
	return newColor;
    }
    return colorTable->indexed_value(it->second).as<Color_O>();
}
#endif

    Color_sp Color_O::create(core::Symbol_sp name, uint rgb)
{
  auto ocolor = gctools::GC<Color_O>::allocate();
    ocolor->setColorName(name);
    ocolor->setColor(rgb);
    return ocolor;
}






void Color_O::initialize()
{
    this->Base::initialize();
    this->_ColorName = nil<core::Symbol_O>();
    this->setRGB(1.0,1.0,1.0);
}




#ifdef XML_ARCHIVE
void	Color_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    if ( node->loading() )
    {
	SIMPLE_ERROR("Never load colors - they are handled by createYourColor");
    }
    node->attributeSymbolIfNotNil("colorName",this->_ColorName);
    node->attribute("color",this->_Color);
}
#endif


void	Color_O::fillDoubleArray(double values[4])
{
    values[0] = this->getRed();
    values[1] = this->getGreen();
    values[2] = this->getBlue();
    values[3] = 0.0;
}

void	Color_O::fillFloatArray(float values[4])
{
    values[0] = this->getRed();
    values[1] = this->getGreen();
    values[2] = this->getBlue();
    values[3] = 0.0;
}

string	Color_O::asString()
{
    stringstream ss;
    ss << "rgb(";
    ss << "\""<<this->_ColorName <<"\",";
    ss << this->getRed() << ",";
    ss << this->getGreen() << ",";
    ss << this->getBlue();
    ss << ")";
    return ss.str();
}


void Color_O::setRGB(double r, double g, double b)
{
uint	ir, ig, ib;
    if ( r < 0.0 || r > 1.0 )
    {
	SIMPLE_ERROR("RGB red value must be between 0.0 and 1.0");
    }
    if ( g < 0.0 || g > 1.0 )
    {
	SIMPLE_ERROR("RGB green value must be between 0.0 and 1.0");
    }
    if ( b < 0.0 || b > 1.0 )
    {
	SIMPLE_ERROR("RGB blue value must be between 0.0 and 1.0");
    }
    ir = 255*r;
    ig = 255*g;
    ib = 255*b;
    LOG("Color bytes r({})g({})b({})" , ir , ig , ib );
    this->_Color = (ir<<16)|(ig<<8)|ib;
    LOG("Put into _Color({})" , this->_Color );
    LOG("Color name({}" , _rep_(this->_ColorName) );
}




    void	Color_O::setColorName(core::Symbol_sp name)
{
  LOG("Setting color name to {}"% _rep_(name) );
    this->_ColorName = name;
}


core::ubyte Color_O::getRedByte()
    {
	core::ubyte ic = ( this->_Color >> 16 )&0xFF;
	return ic;
    }

core::ubyte Color_O::getGreenByte()
    {
	core::ubyte ic = ( this->_Color >> 8 )&0xFF;
	return ic;
    }
core::ubyte Color_O::getBlueByte()
    {
	core::ubyte ic = (this->_Color)&0xFF;
	return ic;
    }

core::ubyte Color_O::getAlphaByte()
    {
	core::ubyte ic = 0xFF;
	return ic;
    }


    double Color_O::getAlpha()
    {
	return 1.0;
    }

double	Color_O::getRed()
{
    uint ic = this->getRedByte();
    double dc = ic/255.0;
    return dc;
}

double	Color_O::getGreen()
{
    uint ic = this->getGreenByte();
    double dc = ic/255.0;
    return dc;
}

double	Color_O::getBlue()
{
    uint ic = this->getBlueByte();
    double dc = ic/255.0;
    return dc;
}



};
