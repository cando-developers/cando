#define	DEBUG_LEVEL_FULL

#include "core/foundation.h"
#include "core/object.h"
#include "core/lisp.h"
#include "core/hashTableEq.h"
#include "symbolTable.h"
#include "core/environment.h"
#include "core/binder.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "core/lisp.h"
#include "color.h"

// last include is wrappers.h
#include "core/wrappers.h"


namespace candoBase
{
    SYMBOL_EXPORT_SC_(CandoBasePkg,colorTable);
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


    void Color_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<Color_O>()
    ;
    Color_O::lisp_initGlobals(_lisp);
}

    void Color_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(CorePkg,Color,"","",_lisp)
    ;
#endif
}


    void Color_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();
	core::HashTableEq_sp colorTable = core::HashTableEq_O::create_default();
	SYMBOL_EXPORT_SC_(CandoBasePkg,colorTable);
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





    Color_O::Color_O(const Color_O& orig) : core::T_O(orig)
{
    this->_ColorName = orig._ColorName;
    this->_Color = orig._Color;
}



    Color_sp Color_O::white(core::Lisp_sp lisp){return Color_O::systemColor(kw::_sym_white);}
    Color_sp Color_O::yellow(core::Lisp_sp lisp){return Color_O::systemColor(kw::_sym_yellow);}
    Color_sp Color_O::red(core::Lisp_sp lisp){return Color_O::systemColor(kw::_sym_red);}
    Color_sp Color_O::green(core::Lisp_sp lisp){return Color_O::systemColor(kw::_sym_green);}
    Color_sp Color_O::blue(core::Lisp_sp lisp){return Color_O::systemColor(kw::_sym_blue);}
    Color_sp Color_O::cyan(core::Lisp_sp lisp){return Color_O::systemColor(kw::_sym_cyan);}
    Color_sp Color_O::magenta(core::Lisp_sp lisp){return Color_O::systemColor(kw::_sym_magenta);}
    Color_sp Color_O::black(core::Lisp_sp lisp){return Color_O::systemColor(kw::_sym_black);}


    Color_sp Color_O::systemColor(core::Symbol_sp sym)
    {
	core::Symbol_sp colorTableSymbol = _sym_colorTable;
	core::Binder_sp colorTable = colorTableSymbol->symbolValue().as<core::Binder_O>();
	return colorTable->lookup(sym).as<Color_O>();
    }


    Color_sp Color_O::getOrDefineSystemColor(core::Symbol_sp sym, uint color)
{
    core::HashTable_sp colorTable = _sym_colorTable->symbolValue().as<core::HashTable_O>();
    Color_sp c = colorTable->gethash(sym);
    if ( c.nilp() ) {
	Color_sp newColor = Color_O::create(sym,color);
	colorTable->setf_gethash(sym,newColor);
	return newColor;
    }
    return c;
}
    
#ifdef XML_ARCHIVE
    Color_sp Color_O::getOrDefineSystemColor(core::ArchiveP node)
{
    core::HashTable_sp colorTable = _sym_colorTable->symbolValue().as<core::HashTable_O>();
    string colorNameString = node->getAttributeString("colorName");
    core::Symbol_sp sym = node->lisp()->internKeyword(colorNameString);
    core::Binder_O::iterator it = colorTable->find(sym);
    if ( it == colorTable->end() )
    {
	uint uicolor = node->getAttributeUInt("color");
	Color_sp newColor = Color_O::create(sym,uicolor,node->lisp());
	colorTable->extend(sym,newColor);
	return newColor;
    }
    return colorTable->indexed_value(it->second).as<Color_O>();
}
#endif

    Color_sp Color_O::create(core::Symbol_sp name, uint rgb)
{
    GC_ALLOCATE(Color_O,ocolor);
    ocolor->setColorName(name);
    ocolor->setColor(rgb);
    return ocolor;
}






void Color_O::initialize()
{
    this->Base::initialize();
    this->_ColorName = _Nil<core::Symbol_O>();
    this->setRGB(1.0,1.0,1.0);
}




#ifdef XML_ARCHIVE
void	Color_O::archiveBase(core::ArchiveP node)
{_G();
    this->Base::archiveBase(node);
    if ( node->loading() )
    {
	SIMPLE_ERROR(BF("Never load colors - they are handled by createYourColor"));
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
{_OF();
uint	ir, ig, ib;
    if ( r < 0.0 || r > 1.0 )
    {
	SIMPLE_ERROR(BF("RGB red value must be between 0.0 and 1.0" ));
    }
    if ( g < 0.0 || g > 1.0 )
    {
	SIMPLE_ERROR(BF("RGB green value must be between 0.0 and 1.0" ));
    }
    if ( b < 0.0 || b > 1.0 )
    {
	SIMPLE_ERROR(BF("RGB blue value must be between 0.0 and 1.0" ));
    }
    ir = 255*r;
    ig = 255*g;
    ib = 255*b;
    LOG(BF("Color bytes r(%d)g(%d)b(%d)") % ir % ig % ib );
    this->_Color = (ir<<16)|(ig<<8)|ib;
    LOG(BF("Put into _Color(%d)") % this->_Color );
    LOG(BF("Color name(%s") % this->_ColorName->__repr__() );
}




    void	Color_O::setColorName(core::Symbol_sp name)
{_G();
    LOG(BF("Setting color name to %s")% name->__repr__() );
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
{_G();
    uint ic = this->getRedByte();
    double dc = ic/255.0;
    return dc;
}

double	Color_O::getGreen()
{_G();
    uint ic = this->getGreenByte();
    double dc = ic/255.0;
    return dc;
}

double	Color_O::getBlue()
{_G();
    uint ic = this->getBlueByte();
    double dc = ic/255.0;
    return dc;
}


EXPOSE_CLASS_AND_GLOBALS(candoBase,Color_O);

};




