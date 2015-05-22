/* *	atomGraph.h * */

#ifndef RingFinder_H  //[
#define RingFinder_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "core/cons.h"
#include "core/bitVector.h"

#include "chemPackage.h"

#include "core/objectSet.fwd.h"// ringFinder.h wants ObjectSet needs objectSet.fwd.h

namespace chem {

SMART(RingFinder);
SMART(AGVertex);
SMART(AGEdge);
SMART(Atom);
SMART(ObjectSet);
SMART(Molecule);
SMART(Matter);

SMART(PathMessage );
class PathMessage_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,PathMessage_O,"PathMessage");
public:
	void initialize();
private:
    RingFinder_wp	_graph;
    core::SimpleBitVector_sp _beep;
    AGVertex_sp	_firstVertex;
    AGEdge_sp	_firstEdge;
    AGVertex_sp	_lastVertex;
public:
    RingFinder_sp getGraph() {_OF(); ANN(this->_graph); return this->_graph.lock(); };
    static PathMessage_sp create(RingFinder_sp graph, AGVertex_sp firstVertex,
    				AGEdge_sp firstEdge, AGVertex_sp lastVertex );

public:
    PathMessage_sp copy();
    void join(PathMessage_sp other);

    AGVertex_sp getFirstVertex();
    AGEdge_sp getFirstEdge();
    AGVertex_sp getLastVertex();
    void updateLastVertex(AGVertex_sp vert);
    void update(AGEdge_sp edge, AGVertex_sp vertex );
    void dump();
    core::Cons_sp getAtoms();
    core::SimpleBitVector_sp	beep();

    PathMessage_O(const PathMessage_O& );

    DEFAULT_CTOR_DTOR(PathMessage_O);
};







SMART(AGVertex );
class AGVertex_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,AGVertex_O,"AGVertex");
private:
    RingFinder_wp		_graph;
    Atom_sp			_atom;
    gctools::Vec0<AGEdge_sp> 		_edges;
    uint			_seenId;
    uint			_backSpan;
    uint			_backCount;
    core::Cons_sp		_sendBuffer;
    core::Cons_sp		_receiveBuffer;
public:
    static AGVertex_sp create(RingFinder_sp graph, Atom_sp atom);
public:

    RingFinder_sp getGraph() {_OF(); ANN(this->_graph);return this->_graph.lock(); };
    void*	getId();
    core::Cons_sp	getConnectedVertices();
    void	dump();
    Atom_sp	getAtom();
    string description() const;
    void addEdge(AGEdge_sp edge);
    void emptySendBuffer();
    void emptyReceiveBuffer();
    void initializeRingSearch();
    void send();
    void acceptMessage(PathMessage_sp msg);
    void receive(uint stage);

    DEFAULT_CTOR_DTOR(AGVertex_O);
};




SMART(AGEdge );
class AGEdge_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,AGEdge_O,"AGEdge");
private:
    RingFinder_wp		_graph;
    uint			_id;
    AGVertex_sp		_vertex1;
    AGVertex_sp		_vertex2;
public:
    static AGEdge_sp create(RingFinder_sp graph, Atom_sp atom1, Atom_sp atom2 );
public:
    RingFinder_sp getGraph() {_OF(); ANN(this->_graph);return this->_graph.lock(); };
    uint getId();
    AGVertex_sp getVertex1();
    AGVertex_sp getVertex2();
    AGVertex_sp otherVertex(AGVertex_O* vert);
    uint getSide(AGVertex_sp vert);
    DEFAULT_CTOR_DTOR(AGEdge_O);
};




SMART(RingFinder );
class RingFinder_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,RingFinder_O,"RingFinder");
private:
    core::HashTable_sp			_vertices;
    gctools::Vec0<AGEdge_sp>		_edges;
    core::HashTableEql_sp               _rings; // m a p<uint,core::Cons_sp>		_rings;
    gctools::Vec0<PathMessage_sp>		_finalRings;
    gctools::Vec0<core::SimpleBitVector_sp>	_gaussian;
public:
    static RingFinder_sp create(Molecule_sp mol);
	/*! Identify the SSSR (Smallest Set of Smallest Rings) for a Matter
	* and set the ring flags for the atoms
	* Return the rings as a Cons of Conses that contain the atoms for each ring
	*/
    static core::Cons_sp identifyRings(Matter_sp matter);
	/*! Identify the SSSR (Smallest Set of Smallest Rings) for a Molecule
	* and set the ring flags for the atoms
	* Return the rings as a Cons of Conses that contain the atoms for each ring
	*/
    static core::Cons_sp identifyRingsInMolecule(Molecule_sp matter);
    /*! Convert the Cons of atoms in one ring into a Cons of bonds that
      join them */
    static core::Cons_sp ringBonds(core::Cons_sp atoms);
      
public:
    void defineForMolecule(Molecule_sp mol);
    AGEdge_sp getEdge(uint i) { return this->_edges[i]; };
    uint getNumberOfEdges();
    uint getNextEdgeId();
//    Atom_sp firstAtom();
    AGVertex_sp vertexForAtom(Atom_sp anAtom);
    void dump();
    core::HashTable_sp& getVertices();
    gctools::Vec0<AGEdge_sp>& getEdges();
    int getNumberOfRingsExpected();
    void addRing(PathMessage_sp ring, uint stage );
    void initializeRingSearch();
    void advanceRingSearch(uint stage);
    void findRings(int numAtoms);
    bool linearlyIndependentRing(PathMessage_sp ring);
    core::Cons_sp getAllRingsAsListsOfAtoms();

    DEFAULT_CTOR_DTOR(RingFinder_O);
};


}; // namespace





TRANSLATE(chem::PathMessage_O);
TRANSLATE(chem::AGVertex_O);
TRANSLATE(chem::AGEdge_O);
TRANSLATE(chem::RingFinder_O);
#endif




