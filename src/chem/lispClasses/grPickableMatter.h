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
#include "addon/vector3.h"
#include "addon/ovector3.h"
#include "core/holder.h"
#include "elements.h"
#include "addon/coordinateArray.h"
#include "bond.h"
#include "addon/highlightObject.h"
#include "grPickableMatter.h"
#include "chemPackage.h"


namespace chem
{


    FORWARD(GrPickableMatter);

    typedef enum { raNone, raPoint, raStick, raBall, raCpk } RenderAtomStyle;
    typedef enum { rbNone, rbLine, rbStick } RenderBondStyle;



    SMART(RenderMatterAtom);
    class RenderMatterAtom_O : public addon::HighlightedObject_O
    {
	LISP_BASE1(addon::HighlightedObject_O);
	LISP_CLASS(chem,ChemPkg,RenderMatterAtom_O,"RenderMatterAtom");
    public:
	
    public:
//	void archiveBase(core::ArchiveP node);
    public:
	core::Symbol_sp	_Style;
	bool		_ShowLabel;
	string		_Label;
	uint 		_NumberOfBonds;
	addon::Color_sp		_Color;
	Atom_sp		_Atom;
    private:
	uint		_BondLineVertexIndex;
    public:
	void	defineForAtom(Atom_sp atom, GrPickableMatter_sp renderMatter);
	Atom_sp getAtom() { return this->_Atom;};
	addon::Color_sp getColor() { return this->_Color;};
	string& getLabel() { return this->_Label;};
	bool getShowLabel() {return this->_ShowLabel;};
	void setBondLineVertexIndex(uint i) { this->_BondLineVertexIndex=i;};
	uint getBondLineVertexIndex() { return this->_BondLineVertexIndex;};

	double distanceToLine(const Vector3& front, const Vector3& back, double& alpha);

	DEFAULT_CTOR_DTOR(RenderMatterAtom_O);
    };

    SMART(RenderMatterBond);
    class RenderMatterBond_O : public addon::HighlightedObject_O
    {
	LISP_BASE1(addon::HighlightedObject_O );
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


    class RenderMatterResidue_O : public addon::HighlightedObject_O
    {
	LISP_BASE1(addon::HighlightedObject_O);
	LISP_CLASS(chem,ChemPkg,RenderMatterResidue_O,"RenderMatterResidue");
    public:
	uint		_PickId;
	Vector3		_Position;
	string		_Label;
	addon::Color_sp	_Color;
    public:
//	void	archiveBase(core::ArchiveP node);
	DEFAULT_CTOR_DTOR(RenderMatterResidue_O);
    };




    SMART(GrPickableMatter);
    class GrPickableMatter_O : public addon::Render_O
    {
	LISP_BASE1(addon::Render_O);
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
	addon::GrPickablePoints_sp		_AtomPoints;
	addon::GrPickableLines_sp		_BondLines;
	addon::GrLabels_sp			_Labels;
    public:
	addon::GrPickablePoints_sp	getAtomPoints() { return this->_AtomPoints;};
	addon::GrPickableLines_sp	getBondLines() { return this->_BondLines;};
	addon::GrLabels_sp		getLabels() { return this->_Labels;};
    public:
	virtual void ownChildren();
	virtual void _walkAndDefineController(addon::RenderController_sp controller, bool suppressSwitch);
	virtual void _walkAndDefineSliderRanges(addon::RenderController_sp controller);

//	virtual void	check();


    private:
//	void	render();
    public:
	static addon::Color_sp elementToColor(Element element, core::Lisp_sp);
    private:
	void generateRenderObjects();
    public:

	void	setFromMatter(Matter_sp obj);
#if 0
	void	setShowLabels(bool t) { this->_ShowLabels = t;};
	bool	getShowLabels() { return this->_ShowLabels;};
#endif
	addon::OVector3_sp centerOfRenderedGeometry();

	virtual addon::Render_sp findChild(core::Cons_sp symbolPath, addon::RenderController_sp controller);


	void accumulateCenterOfGeometry(addon::CenterOfGeometryEnum cogType, addon::RenderController_sp rc, Vector3& sumPos, uint& sumPoints);
	virtual void accumulateBoundingBox(const Matrix& transform, addon::CenterOfGeometryEnum cogType, addon::RenderController_sp rc,
					   addon::BoundingBox_sp bbox);

	void accumulateRenderInformation(stringstream& ss, addon::RenderController_sp cc);


	DEFAULT_CTOR_DTOR(GrPickableMatter_O);
    };





};

TRANSLATE(chem::RenderMatterAtom_O);
TRANSLATE(chem::RenderMatterBond_O);
TRANSLATE(chem::RenderMatterResidue_O);
TRANSLATE(chem::GrPickableMatter_O);
#endif

