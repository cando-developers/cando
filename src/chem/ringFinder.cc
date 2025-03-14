/*
    File: ringFinder.cc
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
#define	DEBUG_LEVEL_NONE

//#define DEBUG_RING_FINDER
//#define DEBUG_RING_MESSAGE


/* This file implements the Smallest Set of Smallest Ring (SSSR) finding algorithm described
   Efficient Exact Solution of the Ring Perception Problem
   Balducci and Pearlman in  J. Chem. Inf. Comput. Sci. Vol 34, No. 4, 1994, p 822-831

   It works well for small molecules <2048 atoms but not for very large systems.
*/

/*
__BEGIN_DOC( candoScript.ringFinder, Ring identification commands and objects)
Commands to identify rings and to manage RingFinder objects.
__END_DOC
*/




#include <clasp/core/foundation.h>
#include <clasp/core/array.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/hashTableEql.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/ringFinder.h>
#include <clasp/core/array.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/adapt/adapters.h>
#include <cando/adapt/objectSet.h>
#include <cando/chem/loop.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/bond.h>
#include <clasp/core/lispDefinitions.h>
#include <clasp/core/wrappers.h>
#include <clasp/core/evaluator.h>

namespace chem {




void PathMessage_O::initialize()
{
    this->Base::initialize();
    this->_beep = nil<core::T_O>();
    this->_firstVertex= nil<AGVertex_O>();
    this->_lastVertex = nil<AGVertex_O>();
    this->_firstEdge = nil<AGEdge_O>();
}


    PathMessage_O::PathMessage_O(const PathMessage_O& orig) : core::CxxObject_O(orig)
{
    this->_graph = orig._graph;
    LOG("About to copy beep" );
    LOG(" beep.nilp() = {}" , orig._beep.nilp()  );
    this->_beep = orig._beep; // shallow copy because we can't allocate in ctor
    LOG("Done copy beep" );
    this->_firstVertex = orig._firstVertex;
    this->_firstEdge = orig._firstEdge;
    this->_lastVertex = orig._lastVertex;
}

PathMessage_sp	PathMessage_O::create(RingFinder_sp graph, AGVertex_sp firstVertex,
    				AGEdge_sp firstEdge, AGVertex_sp lastVertex )
{
  auto  pm  = gctools::GC<PathMessage_O>::allocate_with_default_constructor();
    pm->_graph = graph;
    pm->_beep = core::SimpleBitVector_O::make(graph->getNumberOfEdges());
    pm->_firstVertex = firstVertex;
    pm->_firstEdge = firstEdge;
    pm->_lastVertex = lastVertex;
    pm->_beep->setBit(firstEdge->getId(),1);
    return pm;
}


PathMessage_sp PathMessage_O::copy()
{
  auto  pm  = gctools::GC<PathMessage_O>::copy( *this); // = RP_Copy<PathMessage_O>(this);
    if (this->_beep.nilp()) {
      pm->_beep = nil<core::T_O>();
    } else {
      pm->_beep = core::SimpleBitVector_copy(this->_beep);
    }
    return pm;
}



void PathMessage_O::join(PathMessage_sp other)
{
  core::SimpleBitVector_inPlaceOr(this->_beep,other->_beep);
}

AGVertex_sp PathMessage_O::getFirstVertex()
{
    return this->_firstVertex;
}

AGEdge_sp PathMessage_O::getFirstEdge()
{
    return this->_firstEdge;
}


AGVertex_sp PathMessage_O::getLastVertex()
{
    return this->_lastVertex;
}

void PathMessage_O::updateLastVertex(AGVertex_sp vert)
{
    this->_lastVertex = vert;
}

void PathMessage_O::update(AGEdge_sp edge, AGVertex_sp vertex )
{
    this->_lastVertex = vertex;
    uint ei = edge->getId();
    this->_beep->setBit(ei,1);
}


void PathMessage_O::dump()
{
    Atom_sp a1 = this->_firstVertex->getAtom();
    core::Symbol_sp firstName = a1->getName();
    Atom_sp a2 = this->_lastVertex->getAtom();
    core::Symbol_sp lastName = a2->getName();
    core::clasp_write_string(fmt::format("Ring start: {} end: {}", core::_rep_(firstName) , core::_rep_(lastName)) );
    core::clasp_write_string("Edges: ");
    vector<size_t>::iterator ei;
    vector<size_t> vals;
    core::SimpleBitVector_getOnIndices(this->_beep,vals);
    for ( ei=vals.begin(); ei!=vals.end(); ei++ )
    {
	AGEdge_sp edge = this->getGraph()->getEdge(*ei);
	a1 = edge->getVertex1()->getAtom();
	core::Symbol_sp s1 = a1->getName();
	a2 = edge->getVertex2()->getAtom();
	core::Symbol_sp s2 = a2->getName();
        core::clasp_write_string(fmt::format("    {} - {}"
                                           , s1->symbolName()->get_std_string()
                                           , s2->symbolName()->get_std_string()) );
    }
}


core::List_sp PathMessage_O::getAtoms()
{
    gctools::SmallOrderedSet<Atom_sp>	atoms;
    vector<size_t>::iterator ei;
    vector<size_t> edges;
    SimpleBitVector_getOnIndices(this->_beep,edges);
    LOG("Number of edges = {}" , edges.size()  );
    for ( ei=edges.begin(); ei!=edges.end(); ei++ )
    {
	AGEdge_sp edge = this->getGraph()->getEdge(*ei);
	Atom_sp a1 = edge->getVertex1()->getAtom();
	atoms.insert(a1);
	Atom_sp a2 = edge->getVertex2()->getAtom();
	atoms.insert(a2);
    }
    LOG("Building list" );
    core::List_sp  list = nil<core::T_O>();
    for ( gctools::SmallOrderedSet<Atom_sp>::iterator si=atoms.begin(); si!=atoms.end(); si++ )
    {
	list = core::Cons_O::create(*si,list);
    }
    return list;
}


core::SimpleBitVector_sp	PathMessage_O::beep()
{
    return this->_beep;
}



// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
//
// APVertex
//


AGVertex_sp AGVertex_O::create(RingFinder_sp graph, Atom_sp atom)
{
  auto  v  = gctools::GC<AGVertex_O>::allocate_with_default_constructor();
    v->_graph = graph;
    v->_atom = atom;
    v->_edges.clear();
    v->_seenId = UndefinedUnsignedInt;
    v->emptySendBuffer();
    v->emptyReceiveBuffer();
    return v;
}



core::List_sp	AGVertex_O::getConnectedVertices()
{
    core::List_sp	verts = nil<core::T_O>();
    for ( uint i=0, iEnd(this->_edges.size()); i<iEnd; ++i ) {
	verts = core::Cons_O::create(this->_edges[i],verts);
    }
    return verts;
}


string AGVertex_O::description() const
{
    stringstream ss;
    ss << this->className() << "{";
    ss << this->_atom->description() << "}";
    return ss.str();
}

void	AGVertex_O::dump()
{
  core::clasp_write_string(fmt::format("AGVertex_O {}\n", this->_atom->getName() ));
  core::clasp_write_string(fmt::format("  edges:  "));
  for ( uint i=0, iEnd(this->_edges.size()); i<iEnd; ++i ) {
    core::clasp_write_string(fmt::format("    {} - {} |",
                                         this->_edges[i]->_vertex1->_atom->getName(),
                                         this->_edges[i]->_vertex2->_atom->getName()
                                         ));
  }
  core::clasp_write_string("\n");
}


Atom_sp	AGVertex_O::getAtom()
{
    return this->_atom;
}


void AGVertex_O::addEdge(AGEdge_sp edge)
{
    this->_edges.push_back(edge);
}

void AGVertex_O::emptySendBuffer()
{
    this->_sendBuffer = nil<core::T_O>();
}


void AGVertex_O::emptyReceiveBuffer()
{
    this->_receiveBuffer = nil<core::T_O>();
}


void AGVertex_O::initializeRingSearch()
{
    this->emptySendBuffer();
    for ( uint i=0, iEnd(this->_edges.size()); i<iEnd; ++i ) {
	AGEdge_sp x = this->_edges[i];
	AGVertex_sp nVertex = (x)->otherVertex(this->asSmartPtr());
	PathMessage_sp mess = PathMessage_O::create(this->getGraph(),nVertex,x,nVertex);
	this->_sendBuffer = core::Cons_O::create(mess,this->_sendBuffer);
    }
}


void AGVertex_O::send(int stage)
{
  if ( this->_edges.size() == 0 ) return;
  for ( uint i=0,iEnd(this->_edges.size()); i<iEnd; ++i ) {
    AGEdge_sp edge = this->_edges[i];
    AGVertex_sp neighbor = (edge)->otherVertex(this->asSmartPtr());
#ifdef DEBUG_RING_FINDER
    if (this->_sendBuffer.nilp()) {
      core::clasp_write_string(fmt::format("The send buffer is empty for atom({}) to neighbor({})\n", this->_atom->description(), neighbor->_atom->description() ));
    } else {
      core::clasp_write_string(fmt::format("The send buffer not empty for atom({}) to neighbor({})\n", this->_atom->description(), neighbor->_atom->description() ));
    }
#endif
    for ( auto msgCur : this->_sendBuffer ) {
      PathMessage_sp msg = msgCur->car<PathMessage_O>();
      LOG("Sending message: {}" , _rep_(msg->beep()) );
	    	// If the message most recently came from the neighbor,
		// don't send it back to them
      if ( msg->getLastVertex() == neighbor ) continue;
	    	//
		// If the message originally came from the neighbor 
		// then don't send it back to them
		//
      if ( msg->getFirstVertex() == neighbor ) continue;
      PathMessage_sp newMsg = msg->copy();
      newMsg->update(edge,this->sharedThis<AGVertex_O>());
#ifdef DEBUG_RING_MESSAGE
      if (newMsg->_firstVertex->_atom->getName() == INTERN_(kw,na)) {
        if (
            this->_atom->getName() == INTERN_(kw,cbc) ||
            this->_atom->getName() == INTERN_(kw,cb) ||
            this->_atom->getName() == INTERN_(kw,nb) ||
            this->_atom->getName() == INTERN_(kw,cac) ||
            this->_atom->getName() == INTERN_(kw,ca)
            )  {
          core::clasp_write_string(fmt::format("First :NA : send from {:4} to {:4} -> {}\n", this->_atom->getName(), neighbor->_atom->getName(), _rep_(newMsg->_beep) ));
        }
      }
#endif
      neighbor->acceptMessage(newMsg);
    }
  }
  this->emptySendBuffer();
}

void AGVertex_O::acceptMessage(PathMessage_sp msg)
{
    this->_receiveBuffer = core::Cons_O::create(msg,this->_receiveBuffer);
};


__attribute__((noinline)) PathMessage_sp edgeArrayCopy(PathMessage_sp orig) {
  return orig->copy();
}

__attribute__((noinline)) PathMessage_sp edgeArrayRingCopy(PathMessage_sp orig) {
  return orig->copy();
}

__attribute__((noinline)) void edgeArrayResize(RingFinder_sp graph, gctools::Vec0<PathMessage_sp>& edgeArray) {
  edgeArray.resize(graph->getNumberOfEdges(),nil<PathMessage_O>());
}


void AGVertex_O::receive(uint stage, gctools::Vec0<PathMessage_sp>& edgeArray0, gctools::Vec0<PathMessage_sp>& edgeArray1 )
{
  if ( this->_edges.size() == 0 ) return;
	    //
	    // First cluster all the direct edge messages
	    // and merge them
  RingFinder_sp		graph = this->getGraph();
#if 0
  gctools::Vec0<PathMessage_sp>	edgeArray0, edgeArray1;
  edgeArrayResize(graph,edgeArray0);
  edgeArrayResize(graph,edgeArray1);
#endif
  LOG("Stage({}) Vertex({}) distributing receiveBuffer, there are {} messages" , stage , this->_atom->description().c_str() , this->_receiveBuffer->length() );
  uint addedMessages = 0;
    		// isolated atoms will have empty receive buffers
//  ASSERT(this->_receiveBuffer.notnilp()); //,"The receive buffer is empty for atom("+this->_atom->description()+")!");
#ifdef DEBUG_RING_FINDER
  if (this->_receiveBuffer.nilp()) {
    core::clasp_write_string(fmt::format("The receive buffer is empty for atom({})\n", this->_atom->description() ));
  }
#endif
  if (this->_receiveBuffer.notnilp()) {
    for ( auto msgCur : this->_receiveBuffer ) {
      PathMessage_sp msg = msgCur->car<PathMessage_O>();
      AGEdge_sp edge = msg->getFirstEdge();
      AGVertex_sp vert = msg->getFirstVertex();
#ifdef DEBUG_RING_MESSAGE
      if (msg->_firstVertex->_atom->getName() == INTERN_(kw,na)) {
        if (
            msg->_lastVertex->_atom->getName() == INTERN_(kw,cbc) ||
            msg->_lastVertex->_atom->getName() == INTERN_(kw,cb) ||
            msg->_lastVertex->_atom->getName() == INTERN_(kw,nb) ||
            msg->_lastVertex->_atom->getName() == INTERN_(kw,cac) ||
            msg->_lastVertex->_atom->getName() == INTERN_(kw,ca)
            )  {
        core::clasp_write_string(fmt::format("First :NA : recv from {:4} to {:4} -> {}\n", msg->_lastVertex->_atom->getName(), this->_atom->getName(), _rep_(msg->_beep) ));
        }
      }
#endif
      uint side = edge->getSide(vert);
      uint edgeId = edge->getId();
      if ( side == 0 ) {
        if ( edgeArray0[edgeId].nilp() ) {
          edgeArray0[edgeId] = msg;
          addedMessages++;
        }
      } else {
        if ( edgeArray1[edgeId].nilp() ) {
          edgeArray1[edgeId] = msg;
          addedMessages++;
        }
      }
    }
  }
  LOG("Added {} messages" , addedMessages  );


	    // Now every pair of messages in 
	    //	edgeArray[0][i] and edgeArray[1][i]
	    // share the same first edge but different vertices
	    // They represent odd-membered rings
  uint numEdges;
  { 
    numEdges = graph->getNumberOfEdges();
    for ( uint i=0; i<numEdges; i++ ) {
      LOG("Detecting odd-membered rings i = {}" , i );
		    // detect odd-membered rings
      if ( edgeArray0[i].notnilp() && edgeArray1[i].notnilp() ) {
#ifdef DEBUG_RING_MESSAGE
      if (edgeArray0[i]->_firstVertex->_atom->getName() == INTERN_(kw,na)) {
        core::clasp_write_string(fmt::format("Odd-ring :NA\n"));
      }
#endif
        PathMessage_sp ring = edgeArrayCopy(edgeArray0[i]);
        ring->join(edgeArray1[i]);
        edgeArray0[i] = nil<PathMessage_O>();
        edgeArray1[i] = nil<PathMessage_O>();
        ring->updateLastVertex(this->sharedThis<AGVertex_O>());
        graph->addRing(ring,stage);
      }
    }
  }
	    //
	    // Finally, take all the remaining messages and 
	    // group them by their first vertex
	    // every pair of messages in these groups will be 
	    // vertex collisions
	    //

  core::HashTable_sp vertexDict = core::HashTable_O::createEq();
  {
    LOG("Filling vertexDict" );
    for ( uint i=0; i<numEdges; i++ ) { 
      LOG("Detecting even-membered rings: i = {}" , i );
      if ( edgeArray0[i].notnilp() ) {
        PathMessage_sp msg = edgeArray0[i];
        AGVertex_sp nodeId = msg->getFirstVertex();
        core::List_sp existingMessages = vertexDict->gethash(nodeId,nil<core::T_O>());
        vertexDict->setf_gethash(nodeId,core::Cons_O::create(msg,existingMessages));
      }
      if ( edgeArray1[i].notnilp() ) {
        PathMessage_sp msg = edgeArray1[i];
        AGVertex_sp nodeId = msg->getFirstVertex();
        core::List_sp existingMessages = vertexDict->gethash(nodeId,nil<core::T_O>());
        vertexDict->setf_gethash(nodeId,core::Cons_O::create(msg,existingMessages));
      }
    }
  }
	    //
	    // join all the vertex collision messages in pairs
	    // and add them as rings
	    // messages that aren't collisions are tossed into the
	    // send buffer
	    //
  {
    vertexDict->maphash( [this,&graph,stage] (core::T_sp key, core::T_sp value) {
      core::List_sp list = gc::As<core::List_sp>(value);
		    // If there is more than one PathMessage 
		    // with the same FirstVertex then there is 
		    // a node-collision, combine the messages and
		    // register an even membered ring

      if ( oCdr(list).notnilp() ) {
        for ( auto first : list ) {
          core::Cons_sp second;
          for ( core::List_sp second = oCdr(first); second.notnilp(); second = oCdr(second) ) {
            PathMessage_sp ring = edgeArrayRingCopy(gc::As<PathMessage_sp>(oCar(first)));
            PathMessage_sp secondRing = gc::As<PathMessage_sp>(oCar(second));
            ring->join(secondRing);
            ring->updateLastVertex(this->asSmartPtr());
            graph->addRing(ring,stage);
#ifdef DEBUG_RING_MESSAGE
            if (gc::As<PathMessage_sp>(oCar(first))->_firstVertex->_atom->getName() == INTERN_(kw,na)) {
              core::clasp_write_string(fmt::format("Even-ring :NA -> {} {}\n", ring->_firstVertex->_atom->getName(), _rep_(ring->_beep)));
            }
#endif
          }
        }
      } else {
        this->_sendBuffer = core::Cons_O::create(gc::As<PathMessage_sp>(oCar(list)), this->_sendBuffer);
      }
    } );
  }
  this->emptyReceiveBuffer();
}


AGEdge_sp AGEdge_O::create(RingFinder_sp graph, Atom_sp atom1, Atom_sp atom2 )
{
  auto  edge  = gctools::GC<AGEdge_O>::allocate_with_default_constructor();
    edge->_graph = graph;
    edge->_id = graph->getNextEdgeId();
    if ( atom1->atLowerUniqueAtomOrderThan(atom2) ) {
	edge->_vertex1 = graph->vertexForAtom(atom1);
	edge->_vertex2 = graph->vertexForAtom(atom2);
    } else {
	edge->_vertex2 = graph->vertexForAtom(atom1);
	edge->_vertex1 = graph->vertexForAtom(atom2);
    }
    edge->_vertex1->addEdge(edge);
    edge->_vertex2->addEdge(edge);
    return edge;
}

void AGEdge_O::dump() {
  core::clasp_write_string(fmt::format("AGEdge_O[{}] = {} - {}\n",
                                       this->_id,
                                       this->_vertex1->_atom->getName(),
                                       this->_vertex2->_atom->getName()
                                       ));
}

AGVertex_sp AGEdge_O::getVertex1()
{
    return this->_vertex1;
}

AGVertex_sp AGEdge_O::getVertex2()
{
    return this->_vertex2;
}

uint	AGEdge_O::getId()
{
    return this->_id;
}


// Return the other vertex that this edge points to
AGVertex_sp AGEdge_O::otherVertex(AGVertex_sp firstVertex )
{
    if ( firstVertex == this->_vertex1 )
    {
	return this->_vertex2;
    }
    return this->_vertex1;
}


uint AGEdge_O::getSide(AGVertex_sp vert)
{
    if ( vert == this->_vertex1 )
    {
	return 0;
    }
    return 1;
}


CL_LISPIFY_NAME(make_ring_finder);
CL_DEF_CLASS_METHOD RingFinder_sp RingFinder_O::make(Molecule_sp mol)
{
  auto  graph  = gctools::GC<RingFinder_O>::allocate_with_default_constructor();
    graph->_vertices = core::HashTable_O::create(cl::_sym_eq);
    graph->defineForMolecule(mol);
    return graph;
}


void RingFinder_O::defineForMolecule(Molecule_sp mol)
{
  LOG("Creating RingFinder" );
  Loop loop;
  loop.loopTopGoal(mol,ATOMS);
  this->_vertices->clrhash();
  while ( loop.advanceLoopAndProcess() )
  {
    Atom_sp at = loop.getAtom();
    if ( at.isA<VirtualAtom_O>() ) continue;
    LOG("### Adding atom: {} id:{} to graph" , at->description().c_str() , at.get()  );
//    printf("%s:%d Adding atom as vertex: %s\n",__FILE__,__LINE__,_rep_(at).c_str());
    if ( this->_vertices->gethash(at,nil<T_O>()).notnilp() )
    {
      SIMPLE_ERROR("Non unique atom id");
    }
    AGVertex_sp vert = AGVertex_O::create(this->sharedThis<RingFinder_O>(),at);
    this->_vertices->hash_table_setf_gethash(at,vert);
  }
//  printf("%s:%d Number of vertices: %d\n", __FILE__, __LINE__, this->_vertices->hashTableCount());
  loop.loopTopGoal(mol,BONDS);
  this->_edges.clear();
  while ( loop.advanceLoopAndProcess() )
  {
    Atom_sp a1 = loop.getBondA1();
    Atom_sp a2 = loop.getBondA2();
#if 0
    printf("%s:%d Creating edge between %s - %s\n", __FILE__, __LINE__,
           _rep_(a1).c_str(), _rep_(a2).c_str());
#endif
    AGEdge_sp edge = AGEdge_O::create(this->sharedThis<RingFinder_O>(),a1,a2);
    this->_edges.push_back(edge);
  }
//    printf("%s:%d Number of edges: %lu\n", __FILE__, __LINE__, this->_edges.size());
}

uint RingFinder_O::getNextEdgeId()
{
    return this->_edges.size();
}

#if 0
Atom_sp RingFinder_O::firstAtom()
{
    gctools::SmallMap<Atom_sp,AGVertex_sp>::iterator it;
    it = this->_vertices.begin();
    return it->second->getAtom();
}
#endif

AGVertex_sp RingFinder_O::vertexForAtom(Atom_sp anAtom)
{
    core::T_sp result = this->_vertices->gethash(anAtom,nil<T_O>());
    if (result.nilp()) {
	SIMPLE_ERROR("Could not find atom: {} in atom-vertex map in RingFinder\n" , _rep_(anAtom).c_str() );
    }
    return result.as<AGVertex_O>();
}


class RingFinderVertexMapper : public core::KeyValueMapper
{
public:
    RingFinderVertexMapper() : _Count(0) {};
    int _Count;
    virtual bool mapKeyValue(core::T_sp key, core::T_sp value) {
      gc::As<AGVertex_sp>(value)->dump();
      this->_Count++;
      return true;
    }
};


CL_DEFMETHOD void RingFinder_O::dump()
{
    RingFinderVertexMapper mapper;
    this->_vertices->lowLevelMapHash(&mapper);
    core::clasp_write_string(fmt::format("There are {} atoms\n" , mapper._Count ));
    for ( uint i=0, iEnd(this->_edges.size()); i<iEnd; ++i ) {
      this->_edges[i]->dump();
    }
}


CL_DEFMETHOD core::HashTable_sp RingFinder_O::getVertices()
{
    return this->_vertices;
}

gctools::Vec0<AGEdge_sp>& RingFinder_O::getEdges()
{
    return this->_edges;
}

uint RingFinder_O::getNumberOfEdges()
{
    return this->_edges.size();
}

CL_DEFMETHOD int RingFinder_O::getNumberOfRingsExpected()
{
    LOG("this->_edges.size() = {}" , this->_edges.size()  );
    LOG("this->_vertices->hashTableCount() = {}" , this->_vertices->hashTableCount()  );
    int num = this->_edges.size()-this->_vertices->hashTableCount()+1;
    if (num<0) {
      SIMPLE_ERROR("Predicting {} rings in the molecule - there cannot be a negative number of rings edges: {} and vertices: {}", num, this->_edges.size(), this->_vertices->hashTableCount() );
    }
    LOG("expecting {} rings" , num  );
    return (uint)num;
}



class RingFinderVertexInitializeMapper : public core::KeyValueMapper
{
public:
    virtual bool mapKeyValue(core::T_sp key, core::T_sp value) {
	value.as<AGVertex_O>()->initializeRingSearch();
	return true;
    }
};



void RingFinder_O::initializeRingSearch()
{
    this->_finalRings.clear();
    RingFinderVertexInitializeMapper mapper;
    this->_vertices->lowLevelMapHash(&mapper);
    this->_rings = core::HashTable_O::createEql();
    this->_gaussian.clear();
}


class RingFinderVertexSend : public core::KeyValueMapper
{
  uint _Stage;
public:
  RingFinderVertexSend(uint stage ) : _Stage(stage) {};
public:
    virtual bool mapKeyValue(core::T_sp key, core::T_sp value) {
	value.as<AGVertex_O>()->send(this->_Stage);
#ifdef DEBUG_RING_FINDER
        core::clasp_write_string(fmt::format("Stage: {} sending message from {}\n", this->_Stage, value.as<AGVertex_O>()->_atom->description()));
#endif 
	return true;
    }
};

class RingFinderVertexReceive : public core::KeyValueMapper
{
public:
  uint _Stage;
  size_t _numberOfEdges;
  gctools::Vec0<PathMessage_sp> _edgeArray0;
  gctools::Vec0<PathMessage_sp> _edgeArray1;
  RingFinderVertexReceive(uint stage, size_t numberOfEdges ) : _Stage(stage), _numberOfEdges(numberOfEdges) {
    this->_edgeArray0.resize(numberOfEdges,nil<PathMessage_O>());
    this->_edgeArray1.resize(numberOfEdges,nil<PathMessage_O>());
  };
  virtual bool mapKeyValue(core::T_sp key, core::T_sp value) {
    for ( size_t ii = 0; ii<this->_numberOfEdges; ii++ ) {
      this->_edgeArray0[ii] = nil<PathMessage_O>();
      this->_edgeArray1[ii] = nil<PathMessage_O>();
    }
    value.as<AGVertex_O>()->receive(this->_Stage,this->_edgeArray0,this->_edgeArray1);
    return true;
  }
};


void RingFinder_O::advanceRingSearch(uint stage, size_t numberOfEdges)
{
#ifdef DEBUG_RING_FINDER
  core::clasp_write_string(fmt::format("advanceRingSearch stage: {}\n", stage ));
#endif
#ifdef DEBUG_RING_MESSAGE
  core::clasp_write_string(fmt::format("advanceRingSearch stage: {}\n", stage ));
  core::clasp_write_string("                                                              1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4\n");
  core::clasp_write_string("                                          0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0\n");
#endif
  {
    RingFinderVertexSend sender(stage);
    this->_vertices->lowLevelMapHash(&sender);
  }
  {
    RingFinderVertexReceive receiver(stage, numberOfEdges );
    this->_vertices->lowLevelMapHash(&receiver);
  }
}



CL_LISPIFY_NAME("findRings");
CL_DEFMETHOD void RingFinder_O::findRings(int numAtoms)
{
#if 0
  if (numAtoms >2048) {
    SIMPLE_ERROR("You should not look for rings when there are more than 2048 atoms - the algorithm will blow up memory for {} atoms" , numAtoms );
  }
#endif
  if (numAtoms==0) {
    SIMPLE_ERROR("You tried to find rings in a molecule with zero atoms");
  }
  this->initializeRingSearch();
  //int numberOfRingsExpected = this->getNumberOfRingsExpected();
  //LOG("Number of rings expected = {}" , numberOfRingsExpected  );
  uint stage = 2;
// if we try more than this many times
    // and don't find all the rings then something is wrong
  int trigger = numAtoms*2;
  // Find rings up to 12 atoms
  for (int numSteps = 0; numSteps<6; ++numSteps ) {
//  while ( this->_finalRings.size() < numberOfRingsExpected ) {
    LOG("Looking for rings, stage= {} : trigger[{}] : numAtoms[{}] : numberOfRingsExpected[{}]" , stage , trigger , numAtoms , numberOfRingsExpected  );
    this->advanceRingSearch(stage,this->getNumberOfEdges());
    LOG("After ring search, number of rings found = {}"
        , this->_finalRings.size()  );
    // Stop once we find all the rings we expect to find
    //if (this->_finalRings.size()>=numberOfRingsExpected) break;
    stage++;
  }
}


void RingFinder_O::addRing(PathMessage_sp ring, uint stage)
{
    core::SimpleBitVector_sp beep = ring->beep();
    LOG("Adding ring with beep={}" , _rep_(beep) );
#ifdef DEBUG_RING_MESSAGE
    core::clasp_write_string(fmt::format("addRing with beep: {}\n", _rep_(beep)));
#endif
    core::HashGenerator hg;
    clasp_sxhash(beep,hg);
    core::Fixnum_sp hash = core::clasp_make_fixnum(hg.rawhash());
    core::List_sp ringList = this->_rings->gethash(hash,nil<core::T_O>());
    	// If this new ring is identical to an existing ring then
	// just return
#ifdef DEBUG_RING_MESSAGE
    core::clasp_write_string(fmt::format("comparing rings hash = {}\n", _rep_(hash)));
#endif
    for ( ; ringList.notnilp(); ringList = oCdr(ringList) ) {
      if ( gc::As<PathMessage_sp>(oCar(ringList))->beep()->equal(beep) ) {
#ifdef DEBUG_RING_MESSAGE
        core::clasp_write_string(fmt::format("Found a matching beep - NOT adding ring - returning\n"));
#endif
        return;
      }
    }
    		// Append the ring to the list with this hash
		//
    ringList = this->_rings->gethash(hash,nil<core::T_O>());
    ringList = core::Cons_O::create(ring,ringList);
    this->_rings->setf_gethash(hash,ringList);
    if ( this->linearlyIndependentRing(ring) ) {
#ifdef DEBUG_RING_MESSAGE
      core::clasp_write_string(fmt::format("addRing good ring push_back into ring list\n"));
#endif
	this->_finalRings.push_back(ring);
    } else {
#ifdef DEBUG_RING_MESSAGE
      core::clasp_write_string(fmt::format("Failed linearlyIndependentRing test\n"));
#endif
    }
}


bool RingFinder_O::linearlyIndependentRing(PathMessage_sp ring)
{
#ifdef DEBUG_RING_MESSAGE
  core::clasp_write_string(fmt::format("In linearlyIndependentRing beep             = {}\n", _rep_(ring->beep())));
#endif
  core::SimpleBitVector_sp orig_sbv = ring->beep();
  core::SimpleBitVector_sp beep = core::SimpleBitVector_copy(orig_sbv);
	    //
	    // Swap rows to make the matrix upper triangular
	    //
  uint glast = this->_gaussian.size();
  this->_gaussian.push_back(beep);
//       for z in self._gaussian:
//           print "DEBUG:matrix before: ", z.asString()
//
  if ( this->_gaussian.size() == 1 ) {
#ifdef DEBUG_RING_MESSAGE
    core::clasp_write_string(fmt::format("linearlyIndependentRing returning with gaussian.size==1\n"));
#endif
    return true;
  }

  core::T_sp tleft = core::SimpleBitVector_lowestIndex(this->_gaussian[glast]);
  uint left;
  if (tleft.fixnump()) {
    left = tleft.unsafe_fixnum();
  } else {
    SIMPLE_ERROR("Calculating left, encountered a bitvector full of zeros in the gaussian matrix ");
  }
#ifdef DEBUG_RING_MESSAGE
  core::clasp_write_string(fmt::format("    glast = {}    left = {}\n", glast, left ));
  for ( uint zzz = 0; zzz<this->_gaussian.size(); zzz++ ) {
    core::clasp_write_string(fmt::format("Start linearlyIndependentRing gaussian[{:3}] = {}\n", zzz, _rep_(this->_gaussian[zzz])));
  }
#endif
  for ( uint z = 0; z<glast; z++ ) {
    core::T_sp tgleft = core::SimpleBitVector_lowestIndex(this->_gaussian[z]);
    uint gleft;
    if (tgleft.fixnump()) {
      gleft = tgleft.unsafe_fixnum();
    } else {
      SIMPLE_ERROR("Calculating gleft, encountered a bitvector full of zeros in the gaussian matrix");
    }
#ifdef DEBUG_RING_MESSAGE
    core::clasp_write_string(fmt::format("Checking gaussian[{:3}] = {}\n", z, _rep_(this->_gaussian[z])));
    core::clasp_write_string(fmt::format("    z = {}  gleft = {}\n", z, gleft ));
#endif
    if ( gleft > left ) {
		    // swap rows
	    //print "Swapping rows"
      core::SimpleBitVector_sp temp = this->_gaussian[z];
      this->_gaussian[z] = this->_gaussian[glast];
      this->_gaussian[glast] = temp;
      left = gleft;
#ifdef DEBUG_RING_MESSAGE
      core::clasp_write_string(fmt::format("linearlyIndependentRing swapped {} and {}\n", z, glast ));
#endif
    } else if ( gleft == left ) {
      core::SimpleBitVector_inPlaceXor(this->_gaussian[glast],this->_gaussian[z]);
      tleft = core::SimpleBitVector_lowestIndex(this->_gaussian[glast]);
      if (tleft.fixnump()) {
        left = tleft.unsafe_fixnum();
      } else {
#ifdef DEBUG_RING_MESSAGE
        core::clasp_write_string(fmt::format("After inPlaceXor, linearlyIndependentRing returning false gaussian[{}] = {}\n", glast, _rep_(this->_gaussian[glast])));
#endif
    this->_gaussian.pop_back();
    return false;
      }
#ifdef DEBUG_RING_MESSAGE
      core::clasp_write_string(fmt::format("linearlyIndependentRing nPlaceXor {} and {}\n", glast, z ));
#endif
    }
#ifdef DEBUG_RING_MESSAGE
    for ( uint zzz = 0; zzz<this->_gaussian.size(); zzz++ ) {
      core::clasp_write_string(fmt::format("linearlyIndependentRing gaussian[{}] = {}\n", zzz, _rep_(this->_gaussian[zzz])));
    }
#endif
  }
//  for z in self._gaussian:
//      print "DEBUG:matrix after : ", z.asString()
  if ( core::SimpleBitVector_isZero(this->_gaussian[glast]) ) {
	//print "DEBUG Popped last element"
#ifdef DEBUG_RING_MESSAGE
    core::clasp_write_string(fmt::format("linearlyIndependentRing returning false gaussian[{}] = {}\n", glast, _rep_(this->_gaussian[glast])));
#endif
    this->_gaussian.pop_back();
    return false;
  }
    //print "DEBUG: Returning TRUE"
#ifdef DEBUG_RING_MESSAGE
    core::clasp_write_string(fmt::format("linearlyIndependentRing returning true\n"));
#endif
  return true;
}


CL_LISPIFY_NAME("getAllRingsAsListsOfAtoms");
CL_DEFMETHOD core::List_sp RingFinder_O::getAllRingsAsListsOfAtoms()
{
    core::List_sp lists = nil<core::T_O>();
    gctools::Vec0<PathMessage_sp>::iterator	it;
    uint ridx=0;
    for ( it=this->_finalRings.begin(); it!=this->_finalRings.end(); it++ )
    {
	core::List_sp oneRing = (*it)->getAtoms();
	lists = core::Cons_O::create(oneRing,lists);
	LOG("Ring #{} = size({}) {}" , ridx , core::cl__length(oneRing) , _rep_(oneRing).c_str()  );
	ridx++;
    }
    return lists;
}



core::List_sp RingFinder_O::identifyRingsInMolecule(Molecule_sp molecule)
{
     	//
	// First clear all ring flags
	//
    int numAtoms = 0;
    { 
	Loop loop;
	loop.loopTopGoal(molecule,ATOMS);
	while ( loop.advance() )
	{
	    Atom_sp atom = loop.getAtom();
	    if ( atom->totalBondOrder() > atom->maxTotalBondOrder() )
	    {
		stringstream ss;
		ss << "There is a problem with atom " << atom->description()
		   << " the sum of bond orders is: " << atom->totalBondOrder()
		   << " and it should never be more than " << atom->maxTotalBondOrder();
		SIMPLE_ERROR("{}" , ss.str() );
	    }
	    atom->clearAllRingMembershipFlags();
	    atom->setRingMembershipCount(0);
	    numAtoms++;
	}
    }
    if (numAtoms<3) return nil<core::T_O>();
    RingFinder_sp atomGraph = RingFinder_O::make(molecule);
    {
      atomGraph->findRings(numAtoms);
    }
    core::List_sp rings = nil<core::T_O>();
    {
	rings = atomGraph->getAllRingsAsListsOfAtoms();
	for ( auto curRing : rings ) {
          core::List_sp atoms = oCar(curRing);
          uint ringSize = core::cl__length(atoms);
	    for ( auto atomCons : atoms ) {
              Atom_sp atom = oCar(atomCons).as<Atom_O>();
		atom->setInRingOfSize(ringSize);
		atom->incrementRingMembershipCount();
		LOG("Set {} as part of ring[{}]" , atom->description() , ringSize);
	    }
	}
    }
    return rings;
}






/*
__BEGIN_DOC(candoScript.ringFinder.identifyRings,identifyRings)
\scriptCmd{identifyRings}{matter}

Identify the Smallest Set of Smallest Rings (SSSR) for the Molecule or Aggregate \scriptArg{matter}.
Set the ring membership flags of the atoms that are in rings.
__END_DOC
*/
CL_DEF_CLASS_METHOD core::List_sp RingFinder_O::identifyRings(Matter_sp matter)
{
    if ( matter.isA<Molecule_O>() )
    {
      Molecule_sp mol = (matter).as<Molecule_O>();
	core::List_sp rings = RingFinder_O::identifyRingsInMolecule(mol);
	return rings;
    }
    if ( matter.isA<Aggregate_O>() )
    {
	core::List_sp allRings = nil<core::T_O>();
	Loop molecules;
	molecules.loopTopGoal(matter,MOLECULES);
	while ( molecules.advance() )
	{
	    Molecule_sp mol = molecules.getMolecule();
	    core::List_sp rings = RingFinder_O::identifyRingsInMolecule(mol);
	    // Transfer the rings in rings into allRings
	    while ( rings.notnilp() ) {
		core::List_sp one = rings;
		rings = oCdr(rings);
		one.asCons()->setCdr(allRings);
		allRings = one;
	    }
	}
	return allRings;
    }
    SIMPLE_ERROR("You can only find rings in aggregates or molecules");
}

CL_LISPIFY_NAME(ring-bonds);
DOCGROUP(cando);
CL_DEFUN core::List_sp RingFinder_O::ringBonds(core::List_sp atoms) {
  core::List_sp ringBonds = nil<core::T_O>();
  gctools::SmallOrderedSet<Atom_sp> atomSet;
  {
    
    for (auto cur : atoms) {
      Atom_sp atom = oCar(cur).as<Atom_O>();
      LOG("Atom in ring: {}" , atom->description());
      atomSet.insert(atom);
    }
  }
  int num = 0;
  {
    
    for (auto cur : atoms) {
      Atom_sp atom = oCar(cur).as<Atom_O>();
      VectorBond bonds = atom->getBonds();
      for (gctools::Vec0<Bond_sp>::const_iterator bi = bonds.begin(); bi != bonds.end(); bi++) {
        LOG("Testing {}" , (*bi)->describeOther(atom));
        if (atomSet.count((*bi)->getAtom1()) && atomSet.count((*bi)->getAtom2())) {
          LOG("It's part of the ring!!!");
          core::Cons_sp one = core::Cons_O::create((*bi), ringBonds);
          ringBonds = one;
          num++;
        }
      }
    }
  }
  LOG("There were {} bonds in the ring" , num);
  return ringBonds;
}



};
