/*
    File: atomGrid.cc
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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//

#include "core/common.h"

#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <exception>
#include "core/exceptions.h"
#include "geom/boundingBox.h"
#include "matter.h"
#include "mbbCoreTools.h"
#include "loop.h"
#include "geom/color.h"
#include "geom/render.h"
#include "geom/boundingBox.h"
#include "atomGrid.h"
#include "core/wrappers.h"

#define GRID_PAD 5.0

namespace chem {

void AtomGrid_O::initialize() {
  this->Base::initialize();
  this->xSize = 0;
  this->ySize = 0;
  this->zSize = 0;
}

/*
__BEGIN_DOC(classes.AtomGrid.!class.AtomGrid)
\positional{Matter::matter}
\optionalKeyed{gridResolution}
\optionalKeyed{addRadius}
\optionalKeyed{withinSphere}{List::sphere}
__END_DOC
 */

#if INIT_TO_FACTORIES

#define ARGS_AtomGrid_O_make "(matter grid_resolution add_radius bounding_box)"
#define DECL_AtomGrid_O_make ""
#define DOCS_AtomGrid_O_make "make AtomId args: matter grid_resolution add_radius bounding_box"
AtomGrid_sp AtomGrid_O::make(Matter_sp matter, double gridResolution, double addRadius, geom::BoundingCuboid_sp boundingBox) {
  auto me = gctools::GC<AtomGrid_O>::allocate_with_default_constructor();
  if (boundingBox.notnilp()) {
    me->buildGridWithinBoundingBox(matter, addRadius, gridResolution, boundingBox);
  } else {
    me->buildGrid(matter, addRadius, gridResolution);
  }
  return me;
};

#else

core::T_sp AtomGrid_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp bargs, core::LispPtr lisp) {
  Matter_sp matter = translate::from_object<Matter_sp>::convert(bargs->lookup(ChemPkg, "matter"));
  double stepSize = translate::from_object<core::DoubleFloat_sp>::convert(bargs->lookup(ChemPkg, "gridResolution"))->get();
  double addRadius = translate::from_object<core::DoubleFloat_sp>::convert(bargs->lookup(ChemPkg, "addRadius"))->get();
  geom::BoundingBox_sp bbox = translate::from_object<geom::BoundingBox_sp>::convert(bargs->lookup(ChemPkg, "boundingBox"));
  if (bbox.notnilp()) {
    this->buildGridWithinBoundingBox(matter, addRadius, stepSize, bbox);
  } else {
    this->buildGrid(matter, addRadius, stepSize);
  }
  return _Nil<core::T_O>();
}

#endif

#ifdef XML_ARCHIVE
void AtomGrid_O::archiveBase(core::ArchiveP node) {
  this->Base::archiveBase(node);
  node->attribute("xSize", this->xSize);
  node->attribute("ySize", this->ySize);
  node->attribute("zSize", this->zSize);
  node->attribute("ballRadius", this->ballRadius);
  node->attribute("stepSize", this->stepSize);
  node->attribute("xMin", this->xMin);
  node->attribute("yMin", this->yMin);
  node->attribute("zMin", this->zMin);
  node->attribute("xMax", this->xMax);
  node->attribute("yMax", this->yMax);
  node->attribute("zMax", this->zMax);
  if (node->saving()) {
    LOG("saving grid");
    LOG(" status");
    stringstream sout;
    LOG(" status");
    vector<bool>::iterator it;
    LOG(" status");
    bool done = false;
    LOG(" status");
    it = this->grid.begin();
    LOG(" status done = %d", done);
    while (!done) {
      if (*it) {
        uint oneCount = 0;
        while (it != this->grid.end() && *it) {
          LOG(" status");
          it++;
          LOG(" status");
          oneCount++;
        }
        LOG(" status");
        sout << oneCount << " 1" << std::endl;
        LOG(" entry: %u 1", oneCount);
        if (it == this->grid.end())
          done = true;
        continue;
      }
      uint zeroCount = 0;
      while (it != this->grid.end() && !*it) {
        it++;
        zeroCount++;
      }
      sout << zeroCount << " 0" << std::endl;
      LOG(" entry: %u 1", zeroCount);
      if (it == this->grid.end())
        done = true;
    }
    LOG("Done loop");
    sout << std::endl;
    LOG("Done saving");
    node->setCharacters(sout.str());
  } else {
    LOG("loading grid");
    stringstream sin;
    sin.str(node->characters());
    this->grid.resize(this->xSize * this->ySize * this->zSize, 0);
    vector<bool>::iterator it;
    it = this->grid.begin();
    uint num, val;
    char buffer[1024];
    while (!sin.eof()) {
      sin.getline(buffer, sizeof(buffer));
      string raw = buffer;
      string line = core::trimWhiteSpace(raw);
      LOG("line = |%s|", line.c_str());
      vector<string> parts = core::split(line);
      if (parts.size() == 0)
        continue;
      ASSERT_eq(parts.size(), 2);
      num = atoi(parts[0].c_str());
      val = atoi(parts[1].c_str());
      ASSERT_greaterThan(num, 0);
      LOG("filling %u entries with %u", num, val);
      size_t zzt = (it - this->grid.begin()) + num;
      ASSERT_lt(zzt, this->grid.size() + 1);
      for (uint ui = 0; ui < num; ui++) {
        *it = val;
        it++;
      }
    }
  }
}
#endif

#ifdef XML_ARCHIVE
void AtomGrid_O::archive(core::ArchiveP node) { this->archiveBase(node); }
#endif

void AtomGrid_O::defineMinMaxBox(Matter_sp container, double step, bool useSphere, Vector3 sphereCenter, double sphereRadius) {
  Loop lAtoms;
  double xl, yl, zl;
  double xh, yh, zh;
  Atom_sp a;
  double sphereRadius2;
  Vector3 diff;

  sphereRadius2 = sphereRadius * sphereRadius;
  xl = yl = zl = std::numeric_limits<double>::max();
  xh = yh = zh = -std::numeric_limits<double>::max();
  lAtoms.loopTopGoal(container, ATOMS);
  while (lAtoms.advanceLoopAndProcess()) {
    a = lAtoms.getAtom();
    if (useSphere) {
      diff = a->getPosition() - sphereCenter;
      if (diff.dotProduct(diff) >= sphereRadius2)
        continue;
    }
    xl = MIN(xl, a->getPosition().getX());
    yl = MIN(yl, a->getPosition().getY());
    zl = MIN(zl, a->getPosition().getZ());
    xh = MAX(xh, a->getPosition().getX());
    yh = MAX(yh, a->getPosition().getY());
    zh = MAX(zh, a->getPosition().getZ());
  }
  this->stepSize = step;
  this->xMin = xl - GRID_PAD;
  this->yMin = yl - GRID_PAD;
  this->zMin = zl - GRID_PAD;
  this->xMax = xh + GRID_PAD;
  this->yMax = yh + GRID_PAD;
  this->zMax = zh + GRID_PAD;
  this->xSize = (int)floor((this->xMax - this->xMin) / step) + 1;
  this->ySize = (int)floor((this->yMax - this->yMin) / step) + 1;
  this->zSize = (int)floor((this->zMax - this->zMin) / step) + 1;
  this->grid.resize(this->xSize * this->ySize * this->zSize, false);
}

// Return the index for the position
#define xIndex(ag, xPos) ((int)(floor(((xPos) - ag->xMin) / ag->stepSize)))
#define yIndex(ag, yPos) ((int)(floor(((yPos) - ag->yMin) / ag->stepSize)))
#define zIndex(ag, zPos) ((int)(floor(((zPos) - ag->zMin) / ag->stepSize)))
#define xIndexBounded(ag, xPos) (MIN(MAX(0, xIndex(ag, xPos)), ag->xSize - 1))
#define yIndexBounded(ag, yPos) (MIN(MAX(0, yIndex(ag, yPos)), ag->ySize - 1))
#define zIndexBounded(ag, zPos) (MIN(MAX(0, zIndex(ag, zPos)), ag->zSize - 1))

// Return the index within the grid for the x,y,z indexes
#define gridIndex(gd, xi, yi, zi) (gd->xSize * (gd->ySize * (zi) + (yi)) + (xi))

// Return the coordinate of the index
#define xPos(gd, xi) (gd->xMin + ((xi) * gd->stepSize))
#define yPos(gd, yi) (gd->yMin + ((yi) * gd->stepSize))
#define zPos(gd, zi) (gd->zMin + ((zi) * gd->stepSize))

//
//	paintAtom
//
//	Paint the grid points that are contained within the atom
//
bool AtomGrid_O::paintAtom(const Vector3 &pos, double radius) {
  int xBegin, yBegin, zBegin;
  int xEnd, yEnd, zEnd;
  double xp, yp, zp, rad2, dist2;
  int xi, yi, zi;
  vector<bool>::iterator ip;
  int iPainted;
  int iTotal;
  xBegin = xIndexBounded(this, pos.getX() - radius - 1);
  yBegin = yIndexBounded(this, pos.getY() - radius - 1);
  zBegin = zIndexBounded(this, pos.getZ() - radius - 1);
  xEnd = xIndexBounded(this, pos.getX() + radius + 1);
  yEnd = yIndexBounded(this, pos.getY() + radius + 1);
  zEnd = zIndexBounded(this, pos.getZ() + radius + 1);
  rad2 = radius * radius;
  LOG("Painting atom between: %d,%d,%d and %d,%d,%d", xBegin, yBegin, zBegin, xEnd, yEnd, zEnd);
  LOG("  rad2 = %lf", (rad2));
  iPainted = 0;
  iTotal = 0;
  for (zi = zBegin; zi != zEnd; zi++) {
    for (yi = yBegin; yi != yEnd; yi++) {
      ip = this->grid.begin() + gridIndex(this, xBegin, yi, zi);
      for (xi = xBegin; xi != xEnd; xi++, ip++) {
        xp = xPos(this, xi);
        yp = yPos(this, yi);
        zp = zPos(this, zi);
        xp = xp - pos.getX();
        yp = yp - pos.getY();
        zp = zp - pos.getZ();
        dist2 = xp * xp + yp * yp + zp * zp;
        //
        // If the point is within the radius of the sphere
        // then set it to true
        //
        iTotal++;
        if (dist2 <= rad2) {
          *ip = true;
          iPainted++;
        }
      }
    }
  }
  LOG("  painted %d of %d elements", (iPainted), (iTotal));
  return iPainted > 0;
}

//
//	privateBuildGrid
//
//	Define and fill the grid according to the container.
//	The grid points that are within the vDW radius sphere
//	of each atom are turned on.
//
//	Set the ON_GRID and ATOM_FIXED attributes for atoms that are on the grid
//
void AtomGrid_O::privateBuildGrid(Matter_sp container, double addRadius, double stepSize) {
  Loop lAtoms;
  Atom_sp a;
  Vector3 diff;

  this->grid.resize(this->xSize * this->ySize * this->zSize, false);
  LOG("in build grid");
  this->paint(container, addRadius);
}

void AtomGrid_O::buildGrid(Matter_sp container, double addRadius, double stepSize) {
  Vector3 zero;
  zero.set(0.0, 0.0, 0.0);
  this->ballRadius = addRadius;
  this->defineMinMaxBox(container, stepSize, false, zero, 0.0);
  this->privateBuildGrid(container, addRadius, stepSize);
}

void AtomGrid_O::buildGridWithinSphere(Matter_sp container, double addRadius, double stepSize, Vector3 sphereCenter,
                                       double sphereRadius) {
  this->ballRadius = addRadius;
  this->defineMinMaxBox(container, stepSize, true, sphereCenter, sphereRadius);
  this->privateBuildGrid(container, addRadius, stepSize);
}

void AtomGrid_O::buildGridWithinBoundingBox(Matter_sp container, double addRadius, double stepSize, geom::BoundingBox_sp bbox) {
  ASSERT(bbox->isDefined());
  this->stepSize = stepSize;
  this->xMin = bbox->minCorner().getX();
  this->yMin = bbox->minCorner().getY();
  this->zMin = bbox->minCorner().getZ();
  this->xMax = bbox->maxCorner().getX();
  this->yMax = bbox->maxCorner().getY();
  this->zMax = bbox->maxCorner().getZ();
  this->xSize = (int)floor((this->xMax - this->xMin) / stepSize) + 1;
  this->ySize = (int)floor((this->yMax - this->yMin) / stepSize) + 1;
  this->zSize = (int)floor((this->zMax - this->zMin) / stepSize) + 1;
  this->privateBuildGrid(container, addRadius, stepSize);
}

void AtomGrid_O::paint(Matter_sp matter, double pad) {
  if (matter.notnilp()) {
    if (matter.isA<Atom_O>()) {
      Atom_sp a = matter.as<Atom_O>();
      a->modifyFlags(ATOMFLAG_OFF, ON_GRID | ATOM_FIXED);
      double vdw = vdwRadiusForElement(a->getElement());
      bool onGrid = this->paintAtom(a->getPosition(), vdw + pad);
      if (onGrid) {
        a->modifyFlags(ATOMFLAG_ON, ON_GRID | ATOM_FIXED);
      }
    } else {
      Loop lAtoms;
      lAtoms.loopTopGoal(matter, ATOMS);
      while (lAtoms.advanceLoopAndProcess()) {
        Atom_sp a = lAtoms.getAtom();
        a->modifyFlags(ATOMFLAG_OFF, ON_GRID | ATOM_FIXED);
        double vdw = vdwRadiusForElement(a->getElement());
        bool onGrid = this->paintAtom(a->getPosition(), vdw + pad);
        if (onGrid) {
          a->modifyFlags(ATOMFLAG_ON, ON_GRID | ATOM_FIXED);
        }
      }
    }
  }
}

//
//	dump
//
//	Describe the atomGrid
//
void AtomGrid_O::dump() {
  int iOn, iOff;
  core::clasp_write_string(fmt::format("AtomGrid dump\n"));
  core::clasp_write_string(fmt::format("  size= ({}, {}, {})\n", this->xSize, this->ySize, this->zSize));
  core::clasp_write_string(fmt::format("  min=({}, {}, {})\n", this->xMin, this->yMin, this->zMin));
  core::clasp_write_string(fmt::format("  max=({}, {}, {})\n", this->xMax, this->yMax, this->zMax));
  core::clasp_write_string(fmt::format("  ballRadius={}\n", this->ballRadius));
  core::clasp_write_string(fmt::format("  stepSize={}\n", this->stepSize));
  core::clasp_write_string(fmt::format("  grid Elements={}\n", (int)(this->grid.size())));
  iOn = 0;
  iOff = 0;
  for (uint i = 0; i < this->grid.size(); i++) {
    if (this->grid[i])
      iOn++;
    else
      iOff++;
  }
  core::clasp_write_string(fmt::format("  The grid has {} elements on and {} elements off\n", iOn, iOff));
}

void AtomGrid_O::renderSquare(geom::DisplayList_sp dl, double xn, double yn, double zn, double x1, double y1, double z1, double x2,
                              double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4) {
  Vector3 v1, v2, v3, v4, vn;
  IMPLEMENT_ME();
#if 0
    vn.set(xn,yn,zn);
    v1.set(x1,y1,z1);
    v2.set(x2,y2,z2);
    v3.set(x3,y3,z3);
    v4.set(x4,y4,z4);
    RPGrPolygon polygon = O_GrPolygon::create();
    polygon->setSurfaceNormal(vn);
    polygon->addVertex(v1);
    polygon->addVertex(v2);
    polygon->addVertex(v3);
    polygon->addVertex(v4);
    dl->add(polygon);
#if 0 // if polygons dont work, try lines
    dl->add(geom::GrLine_O::create(v1,v2,1));
    dl->add(geom::GrLine_O::create(v2,v3,1));
    dl->add(geom::GrLine_O::create(v3,v4,1));
    dl->add(geom::GrLine_O::create(v4,v1,1));
#endif
#endif
}

/*!
 * Draw the outline of the grid with points
 */

#ifdef RENDER
geom::Render_sp AtomGrid_O::rendered(core::Cons_sp options) {
  IMPLEMENT_ME();
#if 0
int		xi, yi, zi;
int		gi;
int		gixm, gixp;
int		giym, giyp;
Vector3		v, v1, v2;
double		halfStep;
double		xmin, ymin, zmin;
double		xmax;

    geom::DisplayList_sp	result = geom::DisplayList_O::create();
    iTotal = 0;
    iOn = 0;

    gsize = this->grid.size();
    LOG("Rendering an AtomGrid" );

		//
		// Render the bounding box
		//
    v1.set(this->xMin, this->yMin, this->zMin );
    v2.set(this->xMax, this->yMin, this->zMin );
    result->add(O_GrColor::systemColor(_lisp->symbol(_kw_red)));
    result->add(geom::GrLine_O::create(v1,v2,2));
    v1.set(this->xMin, this->yMin, this->zMin );
    v2.set(this->xMin, this->yMax, this->zMin );
    result->add(geom::GrLine_O::create(v1,v2,2));
    v1.set(this->xMin, this->yMin, this->zMin );
    v2.set(this->xMin, this->yMin, this->zMax );
    result->add(geom::GrLine_O::create(v1,v2,2));

    v1.set(this->xMax, this->yMax, this->zMax );
    v2.set(this->xMax, this->yMax, this->zMin );
    result->add(geom::GrLine_O::create(v1,v2,2));
    v1.set(this->xMax, this->yMax, this->zMax );
    v2.set(this->xMin, this->yMax, this->zMax );
    result->add(geom::GrLine_O::create(v1,v2,2));
    v1.set(this->xMax, this->yMax, this->zMax );
    v2.set(this->xMax, this->yMin, this->zMax );
    result->add(geom::GrLine_O::create(v1,v2,2));

    v1.set(this->xMax, this->yMin, this->zMin );
    v2.set(this->xMax, this->yMax, this->zMin );
    result->add(geom::GrLine_O::create(v1,v2,2));
    v1.set(this->xMax, this->yMin, this->zMin );
    v2.set(this->xMax, this->yMin, this->zMax );
    result->add(geom::GrLine_O::create(v1,v2,2));

    v1.set(this->xMin, this->yMax, this->zMin );
    v2.set(this->xMax, this->yMax, this->zMin );
    result->add(geom::GrLine_O::create(v1,v2,2));
    v1.set(this->xMin, this->yMax, this->zMin );
    v2.set(this->xMin, this->yMax, this->zMax );
    result->add(geom::GrLine_O::create(v1,v2,2));

    v1.set(this->xMin, this->yMin, this->zMax );
    v2.set(this->xMax, this->yMin, this->zMax );
    result->add(geom::GrLine_O::create(v1,v2,2));
    v1.set(this->xMin, this->yMin, this->zMax );
    v2.set(this->xMin, this->yMax, this->zMax );
    result->add(geom::GrLine_O::create(v1,v2,2));

    RPGrColor color = O_GrColor::create(geom::Color_O::systemColor(_lisp->symbol(_kw_magenta)),_lisp);
    result->add(color);
    halfStep = this->stepSize/2.0;
    for ( zi=0; zi!=this->zSize; zi++ ) {
	LOG("Rendering a Z-slice: %d" , (zi ) );
	for ( yi=0; yi!=this->ySize; yi++ ) {
	    for ( xi=0; xi!=this->xSize; xi++ ) {
		iTotal++;
		gi = gridIndex(this,xi,yi,zi);
		if ( this->grid[gi] ) {
		    iOn++;
		    gixm = gridIndex(this,xi-1,yi,zi);
		    gixp = gridIndex(this,xi+1,yi,zi);
		    giym = gridIndex(this,xi,yi-1,zi);
		    giyp = gridIndex(this,xi,yi+1,zi);
		    gizm = gridIndex(this,xi,yi,zi-1);
		    gizp = gridIndex(this,xi,yi,zi+1);
		    xmin = xPos(this,xi)-halfStep;
		    ymin = yPos(this,yi)-halfStep;
		    zmin = zPos(this,zi)-halfStep;
		    xmax = xPos(this,xi)+halfStep;
		    ymax = yPos(this,yi)+halfStep;
		    zmax = zPos(this,zi)+halfStep;
		    edge = (gixm>0 )?(!this->grid[gixm]):true;
		    if ( edge ) {
			this->renderSquare( result,
						-1.0, 0.0, 0.0,
						xmin,ymin,zmin,
						xmin,ymax,zmin,
						xmin,ymax,zmax,
						xmin,ymin,zmax );
		    }
		    edge = (gixp<gsize)?(!this->grid[gixp]):true;
		    if ( edge ) {
			this->renderSquare( result,
						1.0, 0.0, 0.0,
						xmax,ymin,zmin,
						xmax,ymax,zmin,
						xmax,ymax,zmax,
						xmax,ymin,zmax );
		    }
		    edge = (giym>0)?(!this->grid[giym]):true;
		    if ( edge ) {
			this->renderSquare( result,
						0.0, -1.0, 0.0,
						xmin,ymin,zmin,
						xmax,ymin,zmin,
						xmax,ymin,zmax,
						xmin,ymin,zmax );
		    }
		    edge = (giyp<gsize)?(!this->grid[giyp]):true;
		    if ( edge ) {
			this->renderSquare( result,
						0.0, 1.0, 0.0,
						xmin,ymax,zmin,
						xmax,ymax,zmin,
						xmax,ymax,zmax,
						xmin,ymax,zmax );
		    }
		    edge = (gizm>0)?(!this->grid[gizm]):true;
		    if ( edge ) {
			this->renderSquare( result,
						0.0, 0.0, -1.0,
						xmin,ymin,zmin,
						xmax,ymin,zmin,
						xmax,ymax,zmin,
						xmin,ymax,zmin );
		    }
		    edge = (gizp<gsize)?(!this->grid[gizp]):true;
		    if ( edge ) {
			this->renderSquare( result,
						0.0, 0.0, 1.0,
						xmin,ymin,zmax,
						xmax,ymin,zmax,
						xmax,ymax,zmax,
						xmin,ymax,zmax );
		    }
		}
	    }
	}
    }
    LOG( "Total points=%d on=%d\n"% iTotal , iOn );
    return result;
#endif
}
#endif

void AtomGrid_O::invert() {
  for (uint i = 0; i < this->grid.size(); i++) {
    this->grid[i] = !this->grid[i];
  }
}

//
//	collisionAt
//
//	Return true if there will be a collision if we place an
//	atom/residue at the requested position
//
bool AtomGrid_O::collisionAt(const Vector3 &pos) {
  int xi, yi, zi;
  int ii;
  xi = xIndex(this, pos.getX());
  yi = yIndex(this, pos.getY());
  zi = zIndex(this, pos.getZ());
  LOG("collision test at: %d, %d, %d -->", xi, yi, zi);
  if (xi < 0 || xi >= this->xSize)
    goto NOCOLLISION;
  if (yi < 0 || yi >= this->ySize)
    goto NOCOLLISION;
  if (zi < 0 || zi >= this->zSize)
    goto NOCOLLISION;
  LOG("About to get gridIndex");
  ii = gridIndex(this, xi, yi, zi);
  LOG("Got gridIndex=%d" % ii);
  if (this->grid[ii])
    goto COLLISION;
NOCOLLISION:
  LOG("                no collision");
  return false;
COLLISION:
  LOG("                COLLISION");
  return true;
}

}; // namespace chem
