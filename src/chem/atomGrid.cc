/*
    File: atomGrid.h

    Copied from leap/Ambertools source code.

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

#include <algorithm>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/numbers.h>
#include <cando/chem/atom.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atomGrid.h>




namespace chem {


void	AtomGrid_O::fields(core::Record_sp node)
{
  node->/*pod_*/field_if_not_default<double>( INTERN_(kw,originx), this->_OriginX, 0 );
  node->/*pod_*/field_if_not_default<double>( INTERN_(kw,originy), this->_OriginY, 0 );
  node->/*pod_*/field_if_not_default<double>( INTERN_(kw,originz), this->_OriginZ, 0 );
  node->/*pod_*/field_if_not_default<double>( INTERN_(kw,cellSize), this->_CellSize, 0 );
  node->/*pod_*/field_if_not_default<size_t>( INTERN_(kw,nx), this->_Nx, 0 );
  node->/*pod_*/field_if_not_default<size_t>( INTERN_(kw,ny), this->_Ny, 0 );
  node->/*pod_*/field_if_not_default<size_t>( INTERN_(kw,nz), this->_Nz, 0 );
  node->field( INTERN_(kw,cells), this->_Cells);
  node->field( INTERN_(kw,cellSize), this->_CellSize);
  node->field( INTERN_(kw,cellMaxRadius), this->_CellMaxRadius);
  node->field( INTERN_(kw,positions), this->_Positions);
  node->field( INTERN_(kw,effRadii), this->_EffRadii);
  this->Base::fields(node);
}


CL_DEF_CLASS_METHOD AtomGrid_sp AtomGrid_O::makeAtomGrid(NVector_sp positions,
                                     core::SimpleVector_double_sp effRadii,
                                     double cellSize,
                                     double padding )
{
  auto me = gctools::GC<AtomGrid_O>::allocate();
  me->buildAtomGrid( positions, effRadii, cellSize, padding );
  return me;
}


void AtomGrid_O::buildAtomGrid(NVector_sp positions,
                               core::SimpleVector_double_sp effRadii,
                               double cellSize,
                               double padding) {

  const size_t atomCount = effRadii->size();
  if (positions->size() < atomCount * 3) {
    SIMPLE_ERROR("positions size {} < effRadii size * 3 {}", positions->size(), atomCount * 3);
  }

  this->_HeapBounds = core::ComplexVector_double_O::make_vector(
      64, 0.0, core::make_fixnum(0), nil<core::T_O>(), false, core::make_fixnum(0));
  this->_HeapCells = core::ComplexVector_size_t_O::make_vector(
      64, 0, core::make_fixnum(0), nil<core::T_O>(), false, core::make_fixnum(0));
  this->_Positions = positions;
  this->_EffRadii = effRadii;
  
  if (atomCount == 0) {
    this->_Nx = this->_Ny = this->_Nz = 0;
    this->_Cells.clear();
    this->_CellMaxRadius.clear();
    return;
  }

  double minx = std::numeric_limits<double>::infinity();
  double miny = std::numeric_limits<double>::infinity();
  double minz = std::numeric_limits<double>::infinity();
  double maxx = -std::numeric_limits<double>::infinity();
  double maxy = -std::numeric_limits<double>::infinity();
  double maxz = -std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < atomCount; ++i) {
    size_t base = i * 3;
    double x = (*positions)[base + 0];
    double y = (*positions)[base + 1];
    double z = (*positions)[base + 2];
    double r = (*effRadii)[i] + padding;

    minx = std::min(minx, x - r);
    miny = std::min(miny, y - r);
    minz = std::min(minz, z - r);
    maxx = std::max(maxx, x + r);
    maxy = std::max(maxy, y + r);
    maxz = std::max(maxz, z + r);
  }

  this->_CellSize = cellSize;
  this->_OriginX = minx;
  this->_OriginY = miny;
  this->_OriginZ = minz;

  auto spanX = maxx - minx;
  auto spanY = maxy - miny;
  auto spanZ = maxz - minz;

  this->_Nx = std::max<size_t>(1, static_cast<size_t>(std::floor(spanX / cellSize)) + 1);
  this->_Ny = std::max<size_t>(1, static_cast<size_t>(std::floor(spanY / cellSize)) + 1);
  this->_Nz = std::max<size_t>(1, static_cast<size_t>(std::floor(spanZ / cellSize)) + 1);

  size_t cellCount = this->_Nx * this->_Ny * this->_Nz;
  this->_Cells.resize(cellCount, nil<core::T_O>());
  this->_CellMaxRadius.assign(cellCount, 0.0);
  this->_Visited.resize(cellCount,0);
  this->_VisitedFlag = 0;
  auto clampIndex = [](long idx, size_t n) -> size_t {
    if (idx < 0) return 0;
    if (idx >= static_cast<long>(n)) return n - 1;
    return static_cast<size_t>(idx);
  };

  for (size_t i = 0; i < atomCount; ++i) {
    size_t base = i * 3;
    double x = (*positions)[base + 0];
    double y = (*positions)[base + 1];
    double z = (*positions)[base + 2];
    double r = (*effRadii)[i] + padding;

    long ix0 = static_cast<long>(std::floor((x - r - this->_OriginX) / this->_CellSize));
    long ix1 = static_cast<long>(std::floor((x + r - this->_OriginX) / this->_CellSize));
    long iy0 = static_cast<long>(std::floor((y - r - this->_OriginY) / this->_CellSize));
    long iy1 = static_cast<long>(std::floor((y + r - this->_OriginY) / this->_CellSize));
    long iz0 = static_cast<long>(std::floor((z - r - this->_OriginZ) / this->_CellSize));
    long iz1 = static_cast<long>(std::floor((z + r - this->_OriginZ) / this->_CellSize));

    size_t cix0 = clampIndex(ix0, this->_Nx);
    size_t cix1 = clampIndex(ix1, this->_Nx);
    size_t ciy0 = clampIndex(iy0, this->_Ny);
    size_t ciy1 = clampIndex(iy1, this->_Ny);
    size_t ciz0 = clampIndex(iz0, this->_Nz);
    size_t ciz1 = clampIndex(iz1, this->_Nz);

    for (size_t iz = ciz0; iz <= ciz1; ++iz) {
      for (size_t iy = ciy0; iy <= ciy1; ++iy) {
        for (size_t ix = cix0; ix <= cix1; ++ix) {
          size_t idx = this->cellIndex(ix, iy, iz);

          core::T_sp list = this->_Cells[idx];
          // store atom index directly in the list
          core::T_sp node = core::Cons_O::create(core::make_fixnum(static_cast<core::Fixnum>(i)), list);
          this->_Cells[idx] = node;

          if (r > this->_CellMaxRadius[idx]) {
            this->_CellMaxRadius[idx] = r;
          }
        }
      }
    }
  }
#if 0
  core::print(fmt::format("AtomGrid minx,miny,minz = {}, {}, {}\n", minx, miny, minz ));
  core::print(fmt::format("AtomGrid maxx,maxy,maxz = {}, {}, {}\n", maxx, maxy, maxz ));
  core::print(fmt::format("AtomGrid CellSize = {}\n", this->_CellSize ));
  core::print(fmt::format("AtomGrid nx,ny,nz = {}, {}, {}\n", this->_Nx, this->_Ny, this->_Nz ));
#endif
}




size_t AtomGrid_O::cellIndex(size_t ix, size_t iy, size_t iz) const {
  return (ix+(this->_Nx*(iy+(this->_Ny*iz))));
}

void AtomGrid_O::decodeIndex(size_t idx, size_t& cx, size_t& cy, size_t& cz ) const {
  size_t xy = this->_Nx*this->_Ny;
  cz = idx / xy;
  size_t rem = idx - (cz*xy);
  cy = rem/this->_Nx;
  cx = rem - (cy*this->_Nx);
}

double AtomGrid_O::cellMinDist(double x, double y, double z, size_t ix, size_t iy, size_t iz ) const {
  double minx = this->_OriginX+(ix*this->_CellSize);
  double miny = this->_OriginY+(iy*this->_CellSize);
  double minz = this->_OriginZ+(iz*this->_CellSize);
  double maxx = minx + this->_CellSize;
  double maxy = miny + this->_CellSize;
  double maxz = minz + this->_CellSize;
  double dx;
  if (x < minx) dx = minx - x;
  else if ( x > maxx ) dx = x - maxx;
  else dx = 0.0;
  double dy;
  if (y < miny) dy = miny - y;
  else if ( y > maxy ) dy = y - maxy;
  else dy = 0.0;
  double dz;
  if (z < minz) dz = minz - z;
  else if ( z > maxz ) dz = z - maxz;
  else dz = 0.0;
  return sqrt(dx*dx+dy*dy+dz*dz);
}


void heapSwap(core::ComplexVector_double_sp heapBounds, 
              core::ComplexVector_size_t_sp heapCells,
              size_t ii, size_t jj) {
  double temp = (*heapBounds)[ii];
  (*heapBounds)[ii] = (*heapBounds)[jj];
  (*heapBounds)[jj] = temp;
  size_t itemp = (*heapCells)[ii];
  (*heapCells)[ii] = (*heapCells)[jj];
  (*heapCells)[jj] = itemp;
}

void heapPush(core::ComplexVector_double_sp heapBounds,
              core::ComplexVector_size_t_sp heapCells,
              size_t cellIndex,
              double bound) {
  heapBounds->vectorPushExtend(bound);
  heapCells->vectorPushExtend(cellIndex);

  size_t i = heapBounds->fillPointer() - 1;
  while (i > 0) {
    size_t parent = (i - 1) / 2;
    if ((*heapBounds)[i] < (*heapBounds)[parent]) {
      heapSwap(heapBounds, heapCells, i, parent);
      i = parent;
    } else {
      break;
    }
  }
}

bool heapPop(core::ComplexVector_double_sp heapBounds,
             core::ComplexVector_size_t_sp heapCells,
             size_t& outCell,
             double& outBound) {
  size_t n = heapBounds->fillPointer();
  if (n == 0) return false;

  outBound = (*heapBounds)[0];
  outCell = (*heapCells)[0];

  n -= 1;
  if (n > 0) {
    (*heapBounds)[0] = (*heapBounds)[n];
    (*heapCells)[0] = (*heapCells)[n];
  }
  heapBounds->fillPointerSet(n);
  heapCells->fillPointerSet(n);

  size_t i = 0;
  while (true) {
    size_t left = (i * 2) + 1;
    size_t right = left + 1;
    size_t smallest = i;

    if (left < n && (*heapBounds)[left] < (*heapBounds)[smallest]) {
      smallest = left;
    }
    if (right < n && (*heapBounds)[right] < (*heapBounds)[smallest]) {
      smallest = right;
    }
    if (smallest == i) break;

    heapSwap(heapBounds, heapCells, i, smallest);
    i = smallest;
  }

  return true;
}


double AtomGrid_O::sdfGridDouble_(double x, double y, double z, size_t& closestOut) {
  //  CORE::print(fmt::format("Entered sdfGridDouble x,y,z = {}, {}, {}\n", x,y,z));
  if (this->_Cells.size()==0) {
    SIMPLE_ERROR("Cannot call sdfGrid when there are no cells");
  }
  this->_VisitedFlag++; // inc visited flag;
  if (this->_VisitedFlag==0) {
    // Reset visited very 256 calls
    this->_VisitedFlag++;
    std::fill(this->_Visited.begin(), this->_Visited.end(), 0 );
  }

  this->_HeapBounds->fillPointerSet(0);
  this->_HeapCells->fillPointerSet(0);

  auto clampIndex = [](double v, size_t n) -> size_t {
    long idx = static_cast<long>(std::floor(v));
    if (idx < 0) idx = 0;
    if (idx >= static_cast<long>(n)) idx = static_cast<long>(n) - 1;
    return static_cast<size_t>(idx);
  };

  size_t ix = clampIndex((x - this->_OriginX) / this->_CellSize, this->_Nx);
  size_t iy = clampIndex((y - this->_OriginY) / this->_CellSize, this->_Ny);
  size_t iz = clampIndex((z - this->_OriginZ) / this->_CellSize, this->_Nz);

  size_t cells = 0;
  // push-cell equivalent (lambda)
  auto pushCell = [&](size_t& cells, size_t cx, size_t cy, size_t cz) {
    if (cx >= this->_Nx || cy >= this->_Ny || cz >= this->_Nz) return;

    size_t idx = this->cellIndex(cx, cy, cz);
    if (this->_Visited[idx] == this->_VisitedFlag) return;
    this->_Visited[idx] = this->_VisitedFlag;
    double minDist = this->cellMinDist(x, y, z, cx, cy, cz);
    double maxR = this->_CellMaxRadius.empty() ? 0.0 : this->_CellMaxRadius[idx];
    double bound = minDist - maxR;

    //    core::print(fmt::format("pushCell #{} @ {},{},{} bound: {}\n", cells, cx, cy, cz, bound));
    heapPush(this->_HeapBounds, this->_HeapCells, idx, bound);
    cells++;
  };

  // body of %sdf-grid (main loop)
  pushCell( cells, ix, iy, iz);
  double best = 2.4; // std::numeric_limits<double>::infinity();
  size_t bestAtom = SIZE_MAX;
  size_t cellIdx = 0;
  double bound = 0.0;
  while (heapPop(this->_HeapBounds, this->_HeapCells, cellIdx, bound)) {
    if (bound >= best) break;

    core::T_sp cellList = this->_Cells[cellIdx];
    for (core::T_sp cur = cellList; cur.notnilp(); cur = CONS_CDR(cur)) {
      core::T_sp car = CONS_CAR(cur);
      size_t atomIndex = core::clasp_to_size_t(car);

      size_t base = atomIndex * 3;
      if (base+2>=this->_Positions->size()) {
        SIMPLE_ERROR("Out of bounds access of _Positions {} >= {}", base+2, this->_Positions->size());
      }
      double dx = x - (*this->_Positions)[base + 0];
      double dy = y - (*this->_Positions)[base + 1];
      double dz = z - (*this->_Positions)[base + 2];
      double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
      double d = dist - (*this->_EffRadii)[atomIndex];
      if (d < best) {
        best = d;
        bestAtom = atomIndex;
      }
    }

    size_t cx, cy, cz;
    this->decodeIndex(cellIdx, cx, cy, cz);

    if (cx + 1 < this->_Nx) pushCell(cells, cx + 1, cy, cz);
    if (cx > 0)             pushCell(cells, cx - 1, cy, cz);
    if (cy + 1 < this->_Ny) pushCell(cells, cx, cy + 1, cz);
    if (cy > 0)             pushCell(cells, cx, cy - 1, cz);
    if (cz + 1 < this->_Nz) pushCell(cells, cx, cy, cz + 1);
    if (cz > 0)             pushCell(cells, cx, cy, cz - 1);
  }
  //  core::print(fmt::format("sdfGridDouble {},{},{} -> {} cells: {}\n", x,y,z, best, cells));
  closestOut = bestAtom;
  return best;
}

CL_DEFMETHOD
core::T_mv AtomGrid_O::sdfGridDouble(double x, double y, double z) {
  size_t bestAtom;
  double best = this->sdfGridDouble_(x,y,z,bestAtom);
  return Values(core::clasp_make_double_float(best), (bestAtom==SIZE_MAX) ? nil<core::T_O>() : (core::T_sp)core::make_fixnum(bestAtom));
}

CL_DEFMETHOD
core::T_mv AtomGrid_O::sdfGridSingle(float x, float y, float z) {
  size_t bestAtom;
  double best = this->sdfGridDouble_(x,y,z,bestAtom);
  return Values(core::clasp_make_single_float(best), (bestAtom==SIZE_MAX) ? nil<core::T_O>() : (core::T_sp)core::make_fixnum(bestAtom));
}



};

