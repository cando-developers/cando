/*
    File: trajectory.h
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
#ifndef	Trajectory_H //[
#define Trajectory_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/holder.h>
//#include "geom/objectDictionary.fwd.h"
#include <cando/chem/chemPackage.h>
#include <cando/geom/coordinateArray.fwd.h>

namespace chem {

SMART(SimpleVectorCoordinate);
SMART(Matter);
SMART(Atom);

class Trajectory_O;

SMART(TrajectoryFrame);
class TrajectoryFrame_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,TrajectoryFrame_O,"TrajectoryFrame",core::CxxObject_O);
#if INIT_TO_FACTORIES
 public:
    static TrajectoryFrame_sp make();
#else
    DECLARE_INIT();
#endif
friend class Trajectory_O;
public:
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
	geom::SimpleVectorCoordinate_sp	_Coordinates;

public:	// Creation class functions

protected:
    void fillFromMatter(gctools::Vec0<Atom_sp>& atomList );
    void applyToMatter(gctools::Vec0<Atom_sp>& atomList );
public:
	void setCoordinates(geom::SimpleVectorCoordinate_sp ca) { this->_Coordinates = ca;};

	TrajectoryFrame_O( const TrajectoryFrame_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(TrajectoryFrame_O);
};






SMART(Trajectory);
class Trajectory_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,Trajectory_O,"Trajectory",core::CxxObject_O);
#if INIT_TO_FACTORIES
 public:
    static Trajectory_sp make();
#else
    DECLARE_INIT();
#endif
public:
    template <typename T>
        static void expose(T cl)
	{ 
	    cl
	    .def("getMatter",&Trajectory_O::getMatter)
	    .def("numberOfTrajectoryFrames",&Trajectory_O::numberOfTrajectoryFrames)
	    .def("addFrame",&Trajectory_O::addFrame)
	    .def("getTrajectoryFrame",&Trajectory_O::getTrajectoryFrame)
	    .def("applyTrajectoryFrameToMatter",&Trajectory_O::applyTrajectoryFrameToMatter)
	    ;
	}
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
	Matter_sp		_Matter;
    gctools::Vec0<Atom_sp>		_AtomList;
    gctools::Vec0<TrajectoryFrame_sp>	_Frames;
    core::HashTableEq_sp	_Namespace;

private:
	void _setupAtomList(Matter_sp matter);

public:	// Creation class functions

public:
CL_LISPIFY_NAME("getMatter");
CL_DEFMETHOD 	Matter_sp getMatter() const { return this->_Matter;};

	TrajectoryFrame_sp addFrame(Matter_sp matter);

	core::List_sp trajectoryFramesAsList();

	uint numberOfTrajectoryFrames();

	TrajectoryFrame_sp getTrajectoryFrame(uint idx);

	void applyTrajectoryFrameToMatter(TrajectoryFrame_sp f);


#if 0
	bool canRender() { return true; };
#ifdef RENDER
	geom::Render_sp rendered(core::List_sp kargs);
#endif
#endif

	Trajectory_O( const Trajectory_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(Trajectory_O);
};


};
TRANSLATE(chem::TrajectoryFrame_O);
TRANSLATE(chem::Trajectory_O);
#endif //]
