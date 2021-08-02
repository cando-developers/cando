/*
    File: grPickableMatter.h
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


/*
 *
 *	xmlGraphics.h
 *
 *	Generate xmlGraphics strings
 */

#ifndef	GrPickableMatter_H
#define	GrPickableMatter_H


#include <stdio.h>
#include <string>
#include "core/common.h"
#include "candoBase/vector3.h"
#include "candoBase/ovector3.h"
#include "core/holder.h"
#include "elements.h"
#include "candoBase/coordinateArray.h"
#include "bond.h"
#include "candoBase/highlightObject.h"
#include "grPickableMatter.h"
#include "chemPackage.h"


namespace chem
{


    FORWARD(GrPickableMatter);

    typedef enum { raNone, raPoint, raStick, raBall, raCpk } RenderAtomStyle;
    typedef enum { rbNone, rbLine, rbStick } RenderBondStyle;



    SMART(RenderMatterAtom);
    class RenderMatterAtom_O : public candoBase::HighlightedObject_O
    {
	LISP_BASE1(candoBase::HighlightedObject_O);
	LISP_CLASS(chem,ChemPkg,RenderMatterAtom_O,"RenderMatterAtom");
    public:
	
    public:
//	void archiveBase(core::ArchiveP node);
    public:
	core::Symbol_sp	_Style;
	bool		_ShowLabel;
	string		_Label;
	uint 		_NumberOfBonds;
	candoBase::Color_sp		_Color;
	Atom_sp		_Atom;
    private:
	uint		_BondLineVertexIndex;
    public:
	void	defineForAtom(Atom_sp atom, GrPickableMatter_sp renderMatter);
	Atom_sp getAtom() { return this->_Atom;};
	candoBase::Color_sp getColor() { return this->_Color;};
	string& getLabel() { return this->_Label;};
	bool getShowLabel() {return this->_ShowLabel;};
	void setBondLineVertexIndex(uint i) { this->_BondLineVertexIndex=i;};
	uint getBondLineVertexIndex() { return this->_BondLineVertexIndex;};

	double distanceToLine(const Vector3& front, const Vector3& back, double& alpha);

	DEFAULT_CTOR_DTOR(RenderMatterAtom_O);
    };

    SMART(RenderMatterBond);
    class RenderMatterBond_O : public candoBase::HighlightedObject_O
    {
	LISP_BASE1(candoBase::HighlightedObject_O );
	LISP_CLASS(chem,ChemPkg,RenderMatterBond_O,"RenderMatterBond");
    public:
//	void archiveBase(core::ArchiveP node);
    public:
	uint 		_Atom1;
	uint 		_Atom2;
	Vector3		_Center;
	Bond_sp		_Bond;
    public:
	void setup(uint idx1, uint idx2, Atom_sp a1, Atom_sp a2);
	uint getAtom1Index() { return this->_Atom1;};
	uint getAtom2Index() { return this->_Atom2;};
	Vector3& getCenter() { return this->_Center;}
	Bond_sp getBond() { return this->_Bond;};
	DEFAULT_CTOR_DTOR(RenderMatterBond_O);
    };


    class RenderMatterResidue_O : public candoBase::HighlightedObject_O
    {
	LISP_BASE1(candoBase::HighlightedObject_O);
	LISP_CLASS(chem,ChemPkg,RenderMatterResidue_O,"RenderMatterResidue");
    public:
	uint		_PickId;
	Vector3		_Position;
	string		_Label;
	candoBase::Color_sp	_Color;
    public:
//	void	archiveBase(core::ArchiveP node);
	DEFAULT_CTOR_DTOR(RenderMatterResidue_O);
    };




    SMART(GrPickableMatter);
    class GrPickableMatter_O : public candoBase::Render_O
    {
	LISP_BASE1(candoBase::Render_O);
	LISP_CLASS(chem,ChemPkg,GrPickableMatter_O,"GrPickableMatter");
    public:
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
//	bool	loadFinalize(core::ArchiveP node);
    private:
	Matter_sp			_Matter;
        gctools::Vec0<RenderMatterAtom_sp>	_Atoms;
        gctools::Vec0<RenderMatterBond_sp>	_Bonds;
        gctools::Vec0<RenderMatterResidue_sp>	_Residues;
    private:	// do not archive
	candoBase::GrPickablePoints_sp		_AtomPoints;
	candoBase::GrPickableLines_sp		_BondLines;
	candoBase::GrLabels_sp			_Labels;
    public:
	candoBase::GrPickablePoints_sp	getAtomPoints() { return this->_AtomPoints;};
	candoBase::GrPickableLines_sp	getBondLines() { return this->_BondLines;};
	candoBase::GrLabels_sp		getLabels() { return this->_Labels;};
    public:
	virtual void ownChildren();
	virtual void _walkAndDefineController(candoBase::RenderController_sp controller, bool suppressSwitch);
	virtual void _walkAndDefineSliderRanges(candoBase::RenderController_sp controller);

//	virtual void	check();


    private:
//	void	render();
    public:
	static candoBase::Color_sp elementToColor(Element element, core::LispPtr);
    private:
	void generateRenderObjects();
    public:

	void	setFromMatter(Matter_sp obj);
#if 0
	void	setShowLabels(bool t) { this->_ShowLabels = t;};
	bool	getShowLabels() { return this->_ShowLabels;};
#endif
	candoBase::OVector3_sp centerOfRenderedGeometry();

	virtual candoBase::Render_sp findChild(core::Cons_sp symbolPath, candoBase::RenderController_sp controller);


	void accumulateCenterOfGeometry(candoBase::CenterOfGeometryEnum cogType, candoBase::RenderController_sp rc, Vector3& sumPos, uint& sumPoints);
	virtual void accumulateBoundingBox(const Matrix& transform, candoBase::CenterOfGeometryEnum cogType, candoBase::RenderController_sp rc,
					   candoBase::BoundingBox_sp bbox);

	void accumulateRenderInformation(stringstream& ss, candoBase::RenderController_sp cc);


	DEFAULT_CTOR_DTOR(GrPickableMatter_O);
    };





};

TRANSLATE(chem::RenderMatterAtom_O);
TRANSLATE(chem::RenderMatterBond_O);
TRANSLATE(chem::RenderMatterResidue_O);
TRANSLATE(chem::GrPickableMatter_O);
#endif
