/*
    File: virtualSphere.h
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
#ifndef	virtualSphere_H
#define	virtualSphere_H


#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/geom/matrix.h>

#include <cando/chem/chemPackage.h>


namespace chem {

    typedef	enum {vsIdle, vsTranslate, vsRotate, vsScale} VirtualSphereMotionMode;


SMART(VirtualSphere);

class VirtualSphere_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,VirtualSphere_O,"VirtualSphere",core::CxxObject_O);

public:
    static void lisp_initGlobals(core::Lisp_sp lisp);

//    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
	int	_WindowWidth; /*!< Store the width of virtual sphere window */
	int	_WindowHeight; /*!< Store the height of virtual sphere window */
	int	_StartX;	/*!< Store the X coord where mouse down event happened */
	int	_StartY; /*!< Store the Y coord where mouse down event happened */
	int	_MouseX, _MouseY;
	VirtualSphereMotionMode	_MouseMode;	/*!< What is the virtual sphere controlling */
	bool	_Recalculate;

	double	_Scale;		/*!< Store current scale */
	Vector3	_Translate;	/*!< Store translation */
	Vector3	_Center;		/*!< Center position */
	Matrix	_Rotate;		/*!< Rotation only */
	Matrix	_Overall;	/*!< Overall matrix */

public:

public:

	void	setCenter(const Vector3& pos);
	
	void	setTranslate(const Vector3& pos);

	VirtualSphereMotionMode mouseModeFromSymbol(core::Symbol_sp sym);

		/*!
		* Set the x,y coordinate of where the mouse down event
		* started.  Specify the size of the window and what aspect
		* of the view that the VirtualSphere will control.
		*/
	void	setMouseStartSymbol(core::Symbol_sp mode, int x, int y, int w, int h);
	void	setMouseStart( VirtualSphereMotionMode  mode, int x, int y, int w, int h );
	void	setMouseMove( int x, int y );
	void	setMouseXAxisWheelChange(core::Symbol_sp mode, double w);
	void	setMouseYAxisWheelChange(core::Symbol_sp mode, double w);

	void setScaleChange(double yy);


	Matrix	getOverallMatrix();
	double getScale(double minScale, double maxScale);

	void dump();

	void	convertToRelative(int ax, int ay, double& dx, double& dy );

	DEFAULT_CTOR_DTOR(VirtualSphere_O);
};




};



TRANSLATE(chem::VirtualSphere_O);
#endif
