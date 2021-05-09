

// #include "Dimacs.h"
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef DIMACS_H
#define DIMACS_H


#include <utility>
#include <vector>
#include <clasp/core/common.h>
#include <cando/chem/chemPackage.h>

namespace chem {
class Dimacs_O;
};

template <>
struct gctools::GCInfo<chem::Dimacs_O> {
  static bool constexpr CanAllocateWithNoArguments = true;
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = true;
  static GCInfo_policy constexpr Policy = normal;
};


namespace chem {

FORWARD(Dimacs);

    class Dimacs_O : public core::CxxObject_O {
        LISP_CLASS(chem,ChemPkg,Dimacs_O,"Dimacs",core::CxxObject_O);
    public:
        typedef std::pair<int, int> Edge;
    
    protected: public:
      gctools::Vec0<Edge> edges;
      gctools::Vec0<int> degrees;
      size_t numVertices;
      size_t adjacencyMatrixSizeLimit; // TODO: set this from commandline
    
    public:
      static Dimacs_sp make(size_t numVertices, size_t numEdges );
        
    public:
        Dimacs_O(size_t numVertices, size_t numEdges);
      Dimacs_O();
        ~Dimacs_O();

      void add_edge(size_t v1, size_t v2);
      bool load(const char* fname);
      unsigned int getNumVertices() const {return numVertices;}
      unsigned int getNumEdges() const {return edges.size();}
      gctools::Vec0<gctools::Vec0<char> > getAdjacencyMatrix() const;
      gctools::Vec0<int> getDegrees() const {return degrees;}
      void calculateGraphStats(int& maxdegree, int& minDegree, std::vector<float>& degreeHistogram);
      void fixupInternalsForImageSaveLoad(imageSaveLoad::Fixup* fixup) {
        SIMPLE_ERROR(BF("We cannot current image save/load Dimacs's - you could try copying objects in and out of objects that can be image save/loaded"));
      }

    };

};

#endif // DIMACS_H

