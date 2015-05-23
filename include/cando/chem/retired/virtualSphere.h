       
       
//
// (C) 2004 Christian E. Schafmeister
//
#ifndef	virtualSphere_H
#define	virtualSphere_H


#include <clasp/core/common.h>
#include <cando/candoBase/vector3.h>
#include <cando/candoBase/matrix.h>

#include <cando/chem/chemPackage.h>


namespace chem {

    typedef	enum {vsIdle, vsTranslate, vsRotate, vsScale} VirtualSphereMotionMode;


SMART(VirtualSphere);

class VirtualSphere_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,VirtualSphere_O,"VirtualSphere");

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
