/*
    File: atomGrid.h
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
#endif
