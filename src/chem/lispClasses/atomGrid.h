       
       
//
// (C) 2004 Christian E. Schafmeister
//




#ifndef	ATOMGRID_H
#define	ATOMGRID_H

#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <exception>
#include "core/common.h"
#include "matter.h"


#include "chemPackage.h"

#include "candoBase/boundingBox.fwd.h"// atomGrid.h wants BoundingBox needs boundingBox.fwd.h

namespace chem {


typedef	enum { paintConstant, paintExponential } PaintFunction_Enum;



SMART(AtomGrid);
class AtomGrid_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,AtomGrid_O,"AtomGrid");
#if INIT_TO_FACTORIES
 public:
    static AtomGrid_sp make(Matter_sp matter, double gridResolution, double addRadius, candoBase::BoundingBox_sp boundingBox);
#else
    DECLARE_INIT()
#endif
      public:
    //    void	archive(core::ArchiveP node);
    //    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
    int				xSize;
    int				ySize;
    int				zSize;
    double			ballRadius;
    double			stepSize;
    double			xMin;
    double			yMin;
    double			zMin;
    double			xMax;
    double			yMax;
    double			zMax;
    vector<bool>		grid;
public:
private:
void	privateBuildGrid( 	Matter_sp container,
				double	addRadius,
				double	stepSize );
    void	defineMinMaxBox(Matter_sp cont,
				double step,
				bool	useSphere,
				Vector3	sphereCenter,
				double	sphereRadius );
    bool paintAtom( const Vector3& pos, double addRadius );

    void renderSquare( candoBase::DisplayList_sp dl, 
				double xn, double yn, double zn,
				double x1, double y1, double z1,
				double x2, double y2, double z2,
				double x3, double y3, double z3,
				double x4, double y4, double z4 );
public:


    void	buildGrid( Matter_sp container,
			double addRadius,
			double stepSize );

		/*! Build the grid only for atoms within a sphere
		 * This is good for large structures like proteins.
		 */
    void	buildGridWithinSphere( Matter_sp container,
			double addRadius,
			double stepSize,
				       Vector3	sphereCenter,
			double	sphereRadius );

		/*! Build the grid only for atoms within a bounding box
		 * This is good for large structures like proteins.
		 */
    void	buildGridWithinBoundingBox( Matter_sp container,
			double addRadius,
			double stepSize,
			candoBase::BoundingBox_sp bbox );

    void invert();


	/*! Paint the matter into the grid, set every point that
	 * an atom lays on with true.  Pad the vdw radius of each
	 * atom with pad
	 */
    void paint(Matter_sp matter, double pad);
    
    	
    void	dump();

    bool	collisionAt(const Vector3& pos );


    bool canRender() { return true; };
#ifdef RENDER
    candoBase::Render_sp rendered(core::Cons_sp options);
#endif






    DEFAULT_CTOR_DTOR(AtomGrid_O);
};





};
TRANSLATE(chem::AtomGrid_O);
#endif
