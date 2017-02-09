/*
    File: grPickableMatter.cc
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

#include "geom/render.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "core/symbolTable.h"
#include "aggregate.h"
#include "loop.h"
#include "bond.h"
#include "core/array.h"
#include "geom/color.h"
#include "geom/symbolTable.h"
#include "symbolTable.h"
#include "core/cons.h"
//#include "core/xmlLoadArchive.h"
//#include "core/xmlSaveArchive.h"
#include "geom/objectDictionary.h"
#include "geom/ovector3.h"
#include "core/symbolVector.h"
#include "core/environment.h"
#include "geom/integerKeyObjectDictionary.h"
#include "geom/renderController.h"
#include "grPickableMatter.h"
#include "core/wrappers.h"

#include "elements.h"


namespace chem 
{






void RenderMatterAtom_O::defineForAtom(Atom_sp atom,GrPickableMatter_sp renderMatter)
{_OF();
    this->_Atom = atom;
    this->_Style = atom->getPropertyOrDefault(kw::_sym_renderStyle,kw::_sym_ballAndStick).as<core::Symbol_O>();
    this->_ShowLabel = atom->getPropertyOrDefault(kw::_sym_showLabel,_lisp->_false()).isTrue();
    core::Str_sp label = _Nil<atom->getPropertyOrDefault(kw::_sym_label,core::Str_O>()).as<core::Str_O>();
    if ( label.notnilp() )
    {
        this->_Label = label->get();
    } else
    {
        this->_Label = atom->getName();
    }
    this->_NumberOfBonds = atom->numberOfBonds();
    bool colorByElement = atom->getPropertyOrDefault(kw::_sym_colorByElement,_lisp->_true()).isTrue();
    if ( colorByElement )
    {
	this->_Color = (GrPickableMatter_O::elementToColor(atom->getElement(),_lisp));
    } else
    {
	geom::Color_sp color = atom->getPropertyOrDefault(kw::_sym_color,_Nil<geom::Color_O>()).as<geom::Color_O>();
	if ( color.nilp() )
	{
	    this->_Color = geom::Color_O::systemColor(kw::_sym_magenta);
	} else
	{
	    this->_Color = color;
	}
    }
}


#ifdef XML_ARCHIVE
    void RenderMatterAtom_O::archiveBase(core::ArchiveP node)
{_OF();
    this->Base::archiveBase(node);
    node->attribute("style",this->_Style);
    node->attribute("sl",this->_ShowLabel);
    node->attribute("lbl",this->_Label);
    node->attribute("numBonds",this->_NumberOfBonds);
    node->attribute("color",this->_Color);
    node->attribute("atom",this->_Atom);
}
#endif



    /*! Find the distance from this atom to the line defined from p1->p2
      Also return the alpha that gives the closest point on the line by (p1+alpha(p2-p1))
     */

    double RenderMatterAtom_O::distanceToLine(const Vector3& p1, const Vector3& p2, double& alpha )
    {
	Vector3& p3 = this->_Atom->getPositionRef();
	// Solve for alpha
	Vector3 d13 = p1 - p3;
	Vector3	d21 = p2 - p1;
	double num = d13.dotProduct(d21);
	double denom = d21.dotProduct(d21);
	alpha = num/denom;
	// p0 is now the closest point along the line to the point
	Vector3 p0 = p1 + (d21*alpha);
	Vector3 delta = p3 - p0;
	return delta.length();
    }








void RenderMatterBond_O::setup(uint idx1, uint idx2, Atom_sp a1, Atom_sp a2)
{
    this->_Atom1 = idx1;
    this->_Atom2 = idx2;
    this->_Bond = a1->getBondTo(a2);
    Vector3 v1 = a1->getPosition();
    Vector3 v2 = a2->getPosition();
    Vector3 center = (v1.add(v2)).multiplyByScalar(0.5);
    this->_Center = center;
}


#ifdef XML_ARCHIVE
    void RenderMatterBond_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("i1",this->_Atom1);
    node->attribute("i2",this->_Atom2);
    node->archivePlainObject<Vector3>("ctr","V3",this->_Center);
    node->attribute("bond",this->_Bond);
}
#endif












#ifdef XML_ARCHIVE
    void RenderMatterResidue_O::archiveBase(core::ArchiveP node)
    {
	IMPLEMENT_ME();
    }
#endif


geom::Color_sp GrPickableMatter_O::elementToColor(Element element, core::Lisp_sp lisp)
{
    core::Symbol_sp sym = _sym__PLUS_elementColors_PLUS_;
    geom::IntegerKeyObjectDictionary_sp elementColorMap = sym->symbolValue().as<geom::IntegerKeyObjectDictionary_O>();
    if ( elementColorMap->contains(element) )
    {
	return elementColorMap->get(element).as<geom::Color_O>();
    }
    return elementColorMap->get(element_Undefined).as<geom::Color_O>();
}






void	GrPickableMatter_O::initialize()
{_OF();
    this->Base::initialize();
    this->setName(_lisp->internKeyword("matter"));
}


#ifdef XML_ARCHIVE
    void	GrPickableMatter_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("matter",this->_Matter);
    node->archiveVector0("atoms",this->_Atoms);
    node->archiveVector0("bonds",this->_Bonds);
    node->needsFinalization();
}
#endif





#ifdef XML_ARCHIVE
    bool	GrPickableMatter_O::loadFinalize(core::ArchiveP node)
{
    this->Base::loadFinalize(node);
    this->generateRenderObjects();
    return true;
}
#endif


void GrPickableMatter_O::ownChildren()
{
}




#if 0
void GrPickableMatter_O::check()
{
    iterator	it;
    for ( it=this->_PartsDisplayList.begin(); it!=this->_PartsDisplayList.end(); it++ )
    {
	if ( (*it)->isAssignableTo<geom::Render_O>() )
	{
	    LOG(BF( "Child: %s")% (*it)->description().c_str() );
	    if ( (*it)->_InitializationOwner.pointerp() )
	    {
		LOG(BF( "   (*it)->_InitializationOwner.get() = %p")%
		    (*it)->_InitializationOwner.get() );
	    }
	    geom::Render_sp render = (*it).as<geom::Render_O>();
	    LOG(BF( "    (*it)->getParent().get() = %p")% (*it).as<geom::Render_O>()->getParent().get() );
	    if ( render->getParent().get() != this )
	    {
		SIMPLE_ERROR(BF("A child of %s name(%s) specifically(%s) does not have it as a parent") % this->description() % this->getName() % (*it)->description() );
	    }
	    (*it).as<geom::Render_O>()->check();
	}
    }
}
#endif



geom::Render_sp GrPickableMatter_O::findChild(core::Cons_sp symbolPath, geom::RenderController_sp controller)
{_OF();
    if ( !this->isOn(controller) ) return _Nil<geom::Render_O>();
    if ( this->_Name.nilp() )
    {
	IMPLEMENT_ME();
    }
    /* If this is the last node of the path then return ourselves
       if our name matches the path */
    if ( symbolPath->cdr().nilp() )
    {
	if ( symbolPath->car<core::Symbol_O>() == this->_Name ) return this->sharedThis<geom::Render_O>();
    }
    return _Nil<geom::Render_O>();
}






void	GrPickableMatter_O::_walkAndDefineSliderRanges(geom::RenderController_sp controller)
{
}

void	GrPickableMatter_O::_walkAndDefineController(geom::RenderController_sp controller,bool suppressSwitch)
{
    if ( this->getName().notnilp() )
    {
	core::SymbolVector_sp fullPath = this->getFullName();
	controller->createSwitch(fullPath);
    }
}

#if 0
geom::Render_O::iterator GrPickableMatter_O::beginRender(geom::RenderController_sp controller)
{
    core::SymbolVector_sp fullName = this->getFullSwitchName();
    if ( controller->isSwitchOn(fullName) )
    {
	return this->_PartsDisplayList.begin();
    }
    return this->_PartsDisplayList.end();
}

geom::Render_O::iterator GrPickableMatter_O::endRender(geom::RenderController_sp controller)
{
    return this->_PartsDisplayList.end();
}
#endif



geom::OVector3_sp GrPickableMatter_O::centerOfRenderedGeometry()
{
    gctools::Vec0<RenderMatterAtom_sp>::iterator ai;
    Vector3 sum;
    uint num = 0;
    LOG(BF( "Rendering %d atoms")% this->_Atoms.size() );
    for ( ai=this->_Atoms.begin(); ai!=this->_Atoms.end(); ai++ )
    {
	Atom_sp atom = (*ai)->getAtom();
	sum = sum + atom->getPosition();
	num++;
    }
    if ( num == 0 ) return _Nil<geom::OVector3_O>();
    sum = sum.multiplyByScalar(1.0/((double)(num)));
    return geom::OVector3_O::createFromVector3(_lisp,sum);
}


void	GrPickableMatter_O::generateRenderObjects()
{
    this->_AtomPoints = geom::GrPickablePoints_O::create();
    gctools::Vec0<RenderMatterAtom_sp>::iterator ai;
    for ( ai=this->_Atoms.begin(); ai!=this->_Atoms.end(); ai++ )
    {
        LOG(BF( "Atom %s _NumberOfBonds = %d")% (*ai)->_Label % (*ai)->_NumberOfBonds );
	Vector3 pos = (*ai)->getAtom()->getPosition();
	geom::Color_sp color = (*ai)->getColor();
	this->_AtomPoints->appendPickablePoint(pos,color,*ai);
	(*ai)->setBondLineVertexIndex(UndefinedUnsignedInt);
    }
    gctools::Vec0<RenderMatterBond_sp>::iterator	bi;
    this->_BondLines = geom::GrPickableLines_O::create();
    for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
    {
	uint ai1 = (*bi)->getAtom1Index();
	uint ai2 = (*bi)->getAtom2Index();
	RenderMatterAtom_sp ra1 = this->_Atoms[ai1];
	RenderMatterAtom_sp ra2 = this->_Atoms[ai2];
	Vector3 p1 = ra1->getAtom()->getPosition();
	Vector3 p2 = ra2->getAtom()->getPosition();
	geom::Color_sp c1 = ra1->getColor();
	geom::Color_sp c2 = ra2->getColor();
	Vector3& center = (*bi)->getCenter();
	uint bondLineVertexIndex1 = ra1->getBondLineVertexIndex();
	if ( bondLineVertexIndex1 == UndefinedUnsignedInt )
	{
	    bondLineVertexIndex1 = this->_BondLines->appendVertex(ra1->getAtom()->getPosition(),c1);
	    ra1->setBondLineVertexIndex(bondLineVertexIndex1);
	}
	uint bondLineVertexIndex2 = ra2->getBondLineVertexIndex();
	if ( bondLineVertexIndex2 == UndefinedUnsignedInt )
	{
	    bondLineVertexIndex2 = this->_BondLines->appendVertex(ra2->getAtom()->getPosition(),c2);
	    ra2->setBondLineVertexIndex(bondLineVertexIndex2);
	}
	if ( c1 == c2 )
	{
	    this->_BondLines->appendPickableLine(bondLineVertexIndex1,bondLineVertexIndex2,*bi);
	} else
	{
	    uint centerIndex1 = this->_BondLines->appendVertex(center,c1);
	    uint centerIndex2 = this->_BondLines->appendVertex(center,c2);
	    this->_BondLines->appendPickableLine(bondLineVertexIndex1,centerIndex1,*bi);
	    this->_BondLines->appendPickableLine(bondLineVertexIndex2,centerIndex2,*bi);
	}
    }
    {
	this->_Labels = geom::GrLabels_O::create();
	for ( ai=this->_Atoms.begin(); ai!=this->_Atoms.end(); ai++ )
	{
	    if ( (*ai)->getShowLabel() )
	    {
		Vector3& pos = (*ai)->getAtom()->getPositionRef();
		string& label = (*ai)->getLabel();
		geom::Color_sp color = (*ai)->getColor();
		this->_Labels->appendLabel(pos,label,color);
	    }
	}
    }
}


void	GrPickableMatter_O::setFromMatter(Matter_sp matter )
{
    map<Atom_sp,uint>	atomIndices;
    Loop		atomLoop,bondLoop;
    Atom_sp		atom, a1, a2;
    uint		ai1, ai2;
    RenderMatterAtom_sp	grMatterAtom;
    RenderMatterBond_sp	grMatterBond;

    this->_Matter = matter;
    { _BLOCK_TRACE("Storing atoms");
	this->_Atoms.clear();
	atomLoop.loopTopGoal(matter,ATOMS);
	while ( atomLoop.advanceLoopAndProcess() )
	{
	    grMatterAtom = RenderMatterAtom_O::create();
	    atom = atomLoop.getAtom();
	    grMatterAtom->defineForAtom(atom,this->sharedThis<GrPickableMatter_O>());
	    atomIndices[atom] = this->_Atoms.size();
	    this->_Atoms.push_back(grMatterAtom);
	}
    }
    bondLoop.loopTopGoal(matter, BONDS );
    while ( bondLoop.advance() )
    {_BLOCK_TRACE("Storing bonds");
	grMatterBond = RenderMatterBond_O::create();
	a1 = bondLoop.getBondA1();
	a2 = bondLoop.getBondA2();
	if ( atomIndices.count(a1)==1 && atomIndices.count(a2)==1 )
	{
	    ai1 = atomIndices[a1];
	    ai2 = atomIndices[a2];
	    grMatterBond->setup(ai1,ai2,a1,a2);
	    this->_Bonds.push_back(grMatterBond);
	} else
	{
	    SIMPLE_ERROR(BF("A bond was found that contained an atom that hasn't been seen"));
	}
    }
    this->generateRenderObjects();
}

    void GrPickableMatter_O::accumulateCenterOfGeometry(geom::CenterOfGeometryEnum cogType, geom::RenderController_sp rc, Vector3& sumPos, uint& sumPoints )
{
    if ( !this->isOn(rc) ) return;
    this->_AtomPoints->accumulateCenterOfGeometry(cogType,rc,sumPos,sumPoints);
    this->_BondLines->accumulateCenterOfGeometry(cogType,rc,sumPos,sumPoints);
}


    void GrPickableMatter_O::accumulateBoundingBox(const Matrix& transform, geom::CenterOfGeometryEnum cogType, geom::RenderController_sp rc, geom::BoundingBox_sp bbox)
{
    if ( !this->isOn(rc) ) return;
    this->_AtomPoints->accumulateBoundingBox(transform,cogType,rc,bbox);
}

    void GrPickableMatter_O::accumulateRenderInformation(stringstream& ss, geom::RenderController_sp cc)
    {
	// nothing
    }



};
