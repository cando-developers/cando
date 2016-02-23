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
#define	DEBUG_LEVEL_FULL



/*
__BEGIN_DOC( candoScript.ringFinder, Ring identification commands and objects)
Commands to identify rings and to manage RingFinder objects.
__END_DOC
*/




#include <clasp/core/foundation.h>
#include <clasp/core/str.h>
#include <clasp/core/hashTableEql.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/ringFinder.h>
#include <clasp/core/bitVector.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/adapt/adapters.h>
#include <cando/adapt/objectSet.h>
#include <cando/chem/loop.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/bond.h>
#include <clasp/core/vectorObjectsWithFillPtr.h>
#include <clasp/core/lispDefinitions.h>
#include <clasp/core/wrappers.h>


namespace chem {




void PathMessage_O::initialize()
{
    this->Base::initialize();
    this->_beep = _Nil<core::T_O>();
    this->_firstVertex= _Nil<AGVertex_O>();
    this->_lastVertex = _Nil<AGVertex_O>();
    this->_firstEdge = _Nil<AGEdge_O>();
}


    PathMessage_O::PathMessage_O(const PathMessage_O& orig) : core::CxxObject_O(orig)
{_OF();
    this->_graph = orig._graph;
    LOG(BF("About to copy beep") );
    LOG(BF(" beep.nilp() = %d") % orig._beep.nilp()  );
    if ( orig._beep.nilp() ) 
    {
	this->_beep = _Nil<core::T_O>();
	LOG(BF("Setting this->_beep to nil") );
    } else 
    {
	LOG(BF("Copying this->_beep") );
        core::T_sp sbv = orig._beep->deepCopy();
        this->_beep = gc::As<core::SimpleBitVector_sp>(sbv);
    }
    LOG(BF("Done copy beep") );
    this->_firstVertex = orig._firstVertex;
    this->_firstEdge = orig._firstEdge;
    this->_lastVertex = orig._lastVertex;
}

PathMessage_sp	PathMessage_O::create(RingFinder_sp graph, AGVertex_sp firstVertex,
    				AGEdge_sp firstEdge, AGVertex_sp lastVertex )
{
    GC_ALLOCATE(PathMessage_O, pm );
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
    GC_COPY(PathMessage_O, pm , *this); // = RP_Copy<PathMessage_O>(this);
    return pm;
}



void PathMessage_O::join(PathMessage_sp other)
{
    this->_beep->inPlaceOr(other->_beep);
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
    uint firstResidue = a1->getResidueContainedBy()->getId();
    Atom_sp a2 = this->_lastVertex->getAtom();
    core::Symbol_sp lastName = a2->getName();
    uint lastResidue = a2->getResidueContainedBy()->getId();
    _lisp->print(BF("Ring start: %d:%s end: %d:%s") % firstResidue %  firstName
		 % lastResidue % lastName );
    _lisp->print(BF("Edges: "));
    vector<uint>::iterator ei;
    vector<uint> vals;
    this->_beep->getOnIndices(vals);
    for ( ei=vals.begin(); ei!=vals.end(); ei++ )
    {
	AGEdge_sp edge = this->getGraph()->getEdge(*ei);
	a1 = edge->getVertex1()->getAtom();
	uint r1 = a1->getResidueContainedBy()->getId();
	core::Symbol_sp s1 = a1->getName();
	a2 = edge->getVertex2()->getAtom();
	uint r2 = a2->getResidueContainedBy()->getId();
	core::Symbol_sp s2 = a2->getName();
	_lisp->print(BF("    %d:%s - %d:%s") % r1% s1->symbolName()->get() % r2% s2->symbolName()->get() );
    }
}


core::List_sp PathMessage_O::getAtoms()
{
    gctools::SmallOrderedSet<Atom_sp>	atoms;
    vector<uint>::iterator ei;
    vector<uint> edges;
    this->_beep->getOnIndices(edges);
    LOG(BF("Number of edges = %d") % edges.size()  );
    for ( ei=edges.begin(); ei!=edges.end(); ei++ )
    {
	AGEdge_sp edge = this->getGraph()->getEdge(*ei);
	Atom_sp a1 = edge->getVertex1()->getAtom();
	atoms.insert(a1);
	Atom_sp a2 = edge->getVertex2()->getAtom();
	atoms.insert(a2);
    }
    LOG(BF("Building list") );
    core::List_sp  list = _Nil<core::T_O>();
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
    GC_ALLOCATE(AGVertex_O, v );
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
    core::List_sp	verts = _Nil<core::T_O>();
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
    this->_atom->dump();
    for ( auto cur : this->getConnectedVertices() ) {
 	_lisp->prin1(BF("%s ") % cur->car<AGVertex_O>()->getAtom()->getName() );
    }
    _lisp->print(BF(""));
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
    this->_sendBuffer = _Nil<core::T_O>();
}


void AGVertex_O::emptyReceiveBuffer()
{
    this->_receiveBuffer = _Nil<core::T_O>();
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


void AGVertex_O::send()
{
    if ( this->_edges.size() == 0 ) return;
    for ( uint i=0,iEnd(this->_edges.size()); i<iEnd; ++i ) {
	AGEdge_sp edge = this->_edges[i];
	AGVertex_sp neighbor = (edge)->otherVertex(this->asSmartPtr());
	_BLOCK_TRACEF(BF("sending messages to Vertex: %s") % neighbor->description().c_str() );
	for ( auto msgCur : this->_sendBuffer ) {
	    PathMessage_sp msg = msgCur->car<PathMessage_O>();
	    LOG(BF("Sending message: %s") % msg->beep()->__repr__() );
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
	    neighbor->acceptMessage(newMsg);
	}
    }
    this->emptySendBuffer();
}

void AGVertex_O::acceptMessage(PathMessage_sp msg)
{
    this->_receiveBuffer = core::Cons_O::create(msg,this->_receiveBuffer);
};




void AGVertex_O::receive(uint stage)
{
  if ( this->_edges.size() == 0 ) return;
	    //
	    // First cluster all the direct edge messages
	    // and merge them
  gctools::Vec0<PathMessage_sp>	edgeArray0, edgeArray1;
  RingFinder_sp		graph = this->getGraph();
  edgeArray0.resize(graph->getNumberOfEdges(), _Nil<PathMessage_O>());
  edgeArray1.resize(graph->getNumberOfEdges(), _Nil<PathMessage_O>());
  LOG(BF("Stage(%d) Vertex(%s) distributing receiveBuffer, there are %d messages") % stage % this->_atom->description().c_str() % this->_receiveBuffer->length() );
  uint addedMessages = 0;
    		// isolated atoms will have empty receive buffers
  ASSERTP(this->_receiveBuffer.notnilp(),"The receive buffer is empty for atom("+this->_atom->description()+")!");
  for ( auto msgCur : this->_receiveBuffer ) {
    PathMessage_sp msg = msgCur->car<PathMessage_O>();
    _BLOCK_TRACEF(BF("Stage(%d) Vertex(%s) received bitVector: %s") % stage % this->_atom->getName() % msg->beep()->__repr__() );
    AGEdge_sp edge = msg->getFirstEdge();
    AGVertex_sp vert = msg->getFirstVertex();
    uint side = edge->getSide(vert);
    uint edgeId = edge->getId();
    if ( side == 0 ) {
      if ( edgeArray0[edgeId].nilp() )
      {
        edgeArray0[edgeId] = msg;
        addedMessages++;
      }
    } else
    {
      if ( edgeArray1[edgeId].nilp() )
      {
        edgeArray1[edgeId] = msg;
        addedMessages++;
      }
    }
  }
  LOG(BF("Added %d messages") % addedMessages  );


	    // Now every pair of messages in 
	    //	edgeArray[0][i] and edgeArray[1][i]
	    // share the same first edge but different vertices
	    // They represent odd-membered rings
  uint numEdges;
  { _BLOCK_TRACEF(BF("Detecting odd-membered rings numberOfEdges=%d") % graph->getNumberOfEdges() );
    numEdges = graph->getNumberOfEdges();
    for ( uint i=0; i<numEdges; i++ )
    {
      LOG(BF("Detecting odd-membered rings i = %d") % i );
		    // detect odd-membered rings
      if ( edgeArray0[i].notnilp() && edgeArray1[i].notnilp() )
      {
        PathMessage_sp ring = edgeArray0[i]->copy();
        ring->join(edgeArray1[i]);
        edgeArray0[i] = _Nil<PathMessage_O>();
        edgeArray1[i] = _Nil<PathMessage_O>();
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

  core::HashTableEq_sp vertexDict = core::HashTableEq_O::create_default();
  {_BLOCK_TRACE("Detecting even-membered rings");
    LOG(BF("Filling vertexDict") );
    for ( uint i=0; i<numEdges; i++ )
    { 
      LOG(BF("Detecting even-membered rings: i = %d") % i );
      if ( edgeArray0[i].notnilp() ) {
        PathMessage_sp msg = edgeArray0[i];
        AGVertex_sp nodeId = msg->getFirstVertex();
        core::List_sp existingMessages = vertexDict->gethash(nodeId,_Nil<core::T_O>());
        vertexDict->setf_gethash(nodeId,core::Cons_O::create(msg,existingMessages));
      }
      if ( edgeArray1[i].notnilp() ) {
        PathMessage_sp msg = edgeArray1[i];
        AGVertex_sp nodeId = msg->getFirstVertex();
        core::List_sp existingMessages = vertexDict->gethash(nodeId,_Nil<core::T_O>());
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
  {_BLOCK_TRACEF(BF("Detecting even-membered rings"));
    vertexDict->maphash( [this,&graph,stage] (core::T_sp key, core::T_sp value) {
        AGVertex_sp nodeId = gc::As<AGVertex_sp>(key);
        core::List_sp list = gc::As<core::List_sp>(value);
		    // If there is more than one PathMessage 
		    // with the same FirstVertex then there is 
		    // a node-collision, combine the messages and
		    // register an even membered ring

        if ( oCdr(list).notnilp() ) {
          for ( auto first : list ) {
            core::Cons_sp second;
            for ( core::List_sp second = oCdr(first); second.notnilp(); second = oCdr(second) ) {
              PathMessage_sp ring = gc::As<PathMessage_sp>(oCar(first))->copy();
              PathMessage_sp secondRing = gc::As<PathMessage_sp>(oCar(second));
              ring->join(secondRing);
              ring->updateLastVertex(this->asSmartPtr());
              graph->addRing(ring,stage);
            }
          }
        } else
        {
          this->_sendBuffer = core::Cons_O::create(gc::As<PathMessage_sp>(oCar(list)), this->_sendBuffer);
        }
      } );
  }
#if 0 // original code
    map<void*,core::Cons_sp>	vertexDict;
    {_BLOCK_TRACE("Detecting even-membered rings");
	LOG(BF("Filling vertexDict") );
	for ( uint i=0; i<numEdges; i++ )
	{ 
	    LOG(BF("Detecting even-membered rings: i = %d") % i );
	    if ( edgeArray0[i].notnilp() )
	    {
		PathMessage_sp msg = edgeArray0[i];
		void* nodeId = msg->getFirstVertex()->getId();
		if ( vertexDict.count(nodeId) )
		{
		    vertexDict[nodeId] = core::Cons_O::create(msg,vertexDict[nodeId]);
		} else
		{
		    vertexDict[nodeId]=core::Cons_O::create(msg,_Nil<core::T_O>());
		}
	    }
	    if ( edgeArray1[i].notnilp() )
	    {
		PathMessage_sp msg = edgeArray1[i];
		void* nodeId = msg->getFirstVertex()->getId();
		if ( vertexDict.count(nodeId) )
		{
		    vertexDict[nodeId] = core::Cons_O::create(msg,vertexDict[nodeId]);
		} else
		{
		    vertexDict[nodeId]=core::Cons_O::create(msg,_Nil<core::T_O>());
		}
	    }
	}
    }
	    //
	    // join all the vertex collision messages in pairs
	    // and add them as rings
	    // messages that aren't collisions are tossed into the
	    // send buffer
	    //
    {_BLOCK_TRACEF(BF("Detecting even-membered rings"));
	map<void*,core::Cons_sp>::iterator	it;
	for ( it=vertexDict.begin(); it!=vertexDict.end(); it++ )
	{
		    // If there is more than one PathMessage 
		    // with the same FirstVertex then there is 
		    // a node-collision, combine the messages and
		    // register an even membered ring
		    //
	    if ( it->second->cdr().notnilp() )
	    {
		core::Cons_sp first;
		for ( first=it->second; first.notnilp(); first = first->cdr() )
		{
		    core::Cons_sp second;
		    for ( second=first->cdr();second.notnilp();second=second->cdr())
		    {
			PathMessage_sp ring = first->car<PathMessage_O>()->copy();
			ring->join(second->car<PathMessage_O>());
			ring->updateLastVertex(this->sharedThis<AGVertex_O>());
			graph->addRing(ring,stage);
		    }
		}
	    } else
	    {
		this->_sendBuffer = core::Cons_O::create(it->second->car<PathMessage_O>(),
						    this->_sendBuffer,_lisp);
	    }
	}
    }
#endif
  this->emptyReceiveBuffer();
}


  


AGEdge_sp AGEdge_O::create(RingFinder_sp graph, Atom_sp atom1, Atom_sp atom2 )
{
    GC_ALLOCATE(AGEdge_O, edge );
    edge->_graph = graph;
    edge->_id = graph->getNextEdgeId();
    if ( atom1->atLowerAddressThan(atom2) )
    {
	edge->_vertex1 = graph->vertexForAtom(atom1);
	edge->_vertex2 = graph->vertexForAtom(atom2);
    } else
    {
	edge->_vertex2 = graph->vertexForAtom(atom1);
	edge->_vertex1 = graph->vertexForAtom(atom2);
    }
    edge->_vertex1->addEdge(edge);
    edge->_vertex2->addEdge(edge);
    return edge;
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


RingFinder_sp RingFinder_O::create(Molecule_sp mol)
{
    GC_ALLOCATE(RingFinder_O, graph );
    graph->_vertices = core::HashTable_O::create(cl::_sym_eq);
    graph->defineForMolecule(mol);
    return graph;
}


void RingFinder_O::defineForMolecule(Molecule_sp mol)
{
  LOG(BF("Creating RingFinder") );
  Loop loop;
  loop.loopTopGoal(mol,ATOMS);
  this->_vertices->clrhash();
  while ( loop.advanceLoopAndProcess() )
  {
    Atom_sp at = loop.getAtom();
    if ( at.isA<VirtualAtom_O>() ) continue;
    LOG(BF("### Adding atom: %s id:%p to graph") % at->description().c_str() % at.get()  );
//    printf("%s:%d Adding atom as vertex: %s\n",__FILE__,__LINE__,_rep_(at).c_str());
    if ( this->_vertices->gethash(at,_Nil<T_O>()).notnilp() )
    {
      SIMPLE_ERROR(BF("Non unique atom id"));
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
    core::T_sp result = this->_vertices->gethash(anAtom,_Nil<T_O>());
    if (result.nilp()) {
	SIMPLE_ERROR(BF("Could not find atom: %s in atom-vertex map in RingFinder\n") % _rep_(anAtom).c_str() );
    }
    return result.as<AGVertex_O>();
}


class RingFinderVertexMapper : public core::KeyValueMapper
{
public:
    RingFinderVertexMapper() : _Count(0) {};
    int _Count;
    virtual bool mapKeyValue(core::T_sp key, core::T_sp value) {
	value->dump();
	this->_Count++;
	return true;
    }
};


void RingFinder_O::dump()
{
    RingFinderVertexMapper mapper;
    this->_vertices->lowLevelMapHash(&mapper);
    _lisp->print(BF("There are %d atoms") % mapper._Count );
}


core::HashTable_sp& RingFinder_O::getVertices()
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

int RingFinder_O::getNumberOfRingsExpected()
{
    LOG(BF("this->_edges.size() = %d") % this->_edges.size()  );
    LOG(BF("this->_vertices->hashTableCount() = %d") % this->_vertices->hashTableCount()  );
    int num = this->_edges.size()-this->_vertices->hashTableCount()+1;
    ASSERTF(num>=0, BF("There can not be negative number of rings - calculated rings[%d]/edges[%d]/vertices[%d]") % num % this->_edges.size() % this->_vertices->hashTableCount() );
    LOG(BF("expecting %d rings") % num  );
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
    this->_rings = core::HashTableEql_O::create_default();
    this->_gaussian.clear();
}


class RingFinderVertexSend : public core::KeyValueMapper
{
public:
    virtual bool mapKeyValue(core::T_sp key, core::T_sp value) {
	value.as<AGVertex_O>()->send();
	return true;
    }
};

class RingFinderVertexReceive : public core::KeyValueMapper
{
public:
    uint _Stage;
    RingFinderVertexReceive(uint stage) : _Stage(stage) {};
    virtual bool mapKeyValue(core::T_sp key, core::T_sp value) {
	value.as<AGVertex_O>()->receive(this->_Stage);
	return true;
    }
};


void RingFinder_O::advanceRingSearch(uint stage)
{
    {_BLOCK_TRACE("send stage");
	RingFinderVertexSend sender;
	this->_vertices->lowLevelMapHash(&sender);
    }
    {_BLOCK_TRACE("receive stage");
	RingFinderVertexReceive receiver(stage);
	this->_vertices->lowLevelMapHash(&receiver);
    }
}



CL_LISPIFY_NAME("findRings");
CL_DEFMETHOD void RingFinder_O::findRings(int numAtoms)
{
    ASSERTF(numAtoms>0,BF("You tried to find rings in a molecule with zero atoms"));
    this->initializeRingSearch();
    int numberOfRingsExpected = this->getNumberOfRingsExpected();
    LOG(BF("Number of rings expected = %d") % numberOfRingsExpected  );
    uint stage = 2;
// if we try more than this many times
    // and don't find all the rings then something is wrong
    int trigger = numAtoms*2;
    while ( this->_finalRings.size() < numberOfRingsExpected )
    {
	LOG(BF("Looking for rings, stage= %d : trigger[%d] : numAtoms[%d] : numberOfRingsExpected[%d]") % stage % trigger % numAtoms % numberOfRingsExpected  );
        this->advanceRingSearch(stage);
	LOG(BF("After ring search, number of rings found = %d")
	    % this->_finalRings.size()  );
	stage++;
	trigger--;
	if ( trigger<=0 )
	{
	    SIMPLE_ERROR(BF("We advanced the ring search way beyond the number of times we should have needed to - there are %d atoms and %d rings expected and we advanced the search %d times") % numAtoms % numberOfRingsExpected % stage );
	}
    }
}


void RingFinder_O::addRing(PathMessage_sp ring, uint stage)
{
    core::SimpleBitVector_sp beep = ring->beep();
    LOG(BF("Adding ring with beep=%s") % beep->__repr__() );
    core::HashGenerator hg;
    clasp_sxhash(beep,hg);
    core::Fixnum_sp hash = core::clasp_make_fixnum(hg.hash());
    core::List_sp ringList = this->_rings->gethash(hash,_Nil<core::T_O>());
    	// If this new ring is identical to an existing ring then
	// just return
    for ( ; ringList.notnilp(); ringList = oCdr(ringList) )
    {
      if ( gc::As<PathMessage_sp>(oCar(ringList))->beep()->equal(beep) ) return;
    }
    		// Append the ring to the list with this hash
		//
    ringList = this->_rings->gethash(hash,_Nil<core::T_O>());
    ringList = core::Cons_O::create(ring,ringList);
    this->_rings->setf_gethash(hash,ringList);
    if ( this->linearlyIndependentRing(ring) ) {
	this->_finalRings.push_back(ring);
    }
}


bool RingFinder_O::linearlyIndependentRing(PathMessage_sp ring)
{
    core::SimpleBitVector_sp beep = ring->beep()->deepCopy().as<core::SimpleBitVector_O>();
	    //
	    // Swap rows to make the matrix upper triangular
	    //
    uint glast = this->_gaussian.size();
    this->_gaussian.push_back(beep);
//       for z in self._gaussian:
//           print "DEBUG:matrix before: ", z.asString()
//
    if ( this->_gaussian.size() == 1 ) return true;

	    //
	    // figure where to insert the 
	    //
    uint left = this->_gaussian[glast]->lowestIndex();
    for ( uint z = 0; z<glast; z++ )
    {
	uint gleft = this->_gaussian[z]->lowestIndex();
	if ( gleft > left )
	{
		    // swap rows
	    //print "Swapping rows"
	    core::SimpleBitVector_sp temp = this->_gaussian[z];
	    this->_gaussian[z] = this->_gaussian[glast];
	    this->_gaussian[glast] = temp;
	    left = gleft;
	} else if ( gleft == left )
	{
	    this->_gaussian[glast]->inPlaceXor(this->_gaussian[z]);
	    left = this->_gaussian[glast]->lowestIndex();
	}
    }
//  for z in self._gaussian:
//      print "DEBUG:matrix after : ", z.asString()
    if ( this->_gaussian[glast]->isZero() )
    {
	this->_gaussian.pop_back();
	//print "DEBUG Popped last element"
	return false;
    }
    //print "DEBUG: Returning TRUE"
    return true;
}


CL_LISPIFY_NAME("getAllRingsAsListsOfAtoms");
CL_DEFMETHOD core::List_sp RingFinder_O::getAllRingsAsListsOfAtoms()
{
    core::List_sp lists = _Nil<core::T_O>();
    gctools::Vec0<PathMessage_sp>::iterator	it;
    uint ridx=0;
    for ( it=this->_finalRings.begin(); it!=this->_finalRings.end(); it++ )
    {
	core::List_sp oneRing = (*it)->getAtoms();
	lists = core::Cons_O::create(oneRing,lists);
	LOG(BF("Ring #%d = size(%d) %s") % ridx % core::cl__length(oneRing) % oneRing->__repr__().c_str()  );
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
    { _BLOCK_TRACE("Clearing ring membership");
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
		SIMPLE_ERROR(BF("%s")%ss.str() );
	    }
	    atom->clearAllRingMembershipFlags();
	    atom->setRingMembershipCount(0);
	    numAtoms++;
	}
    }
    RingFinder_sp atomGraph = RingFinder_O::create(molecule);
    {_BLOCK_TRACE("Looking for rings");
	atomGraph->findRings(numAtoms);
    }
    core::List_sp rings = _Nil<core::T_O>();
    {_BLOCK_TRACE("Assigning ring membership");
	rings = atomGraph->getAllRingsAsListsOfAtoms();
	for ( auto curRing : rings ) {
          core::List_sp atoms = oCar(curRing);
          uint ringSize = core::cl__length(atoms);
	    for ( auto atomCons : atoms ) {
              Atom_sp atom = oCar(atomCons).as<Atom_O>();
		atom->setInRingOfSize(ringSize);
		atom->incrementRingMembershipCount();
		LOG(BF("Set %s as part of ring[%d]") % atom->description() % ringSize);
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
CL_LISPIFY_NAME(identify-rings);
CL_DEFUN core::List_sp RingFinder_O::identifyRings(Matter_sp matter)
{
    if ( matter.isA<Molecule_O>() )
    {
        Molecule_sp mol = downcast<Molecule_O>(matter);
	core::List_sp rings = RingFinder_O::identifyRingsInMolecule(mol);
	return rings;
    }
    if ( matter.isA<Aggregate_O>() )
    {
	core::List_sp allRings = _Nil<core::T_O>();
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
    SIMPLE_ERROR(BF("You can only find rings in aggregates or molecules"));
}

CL_LISPIFY_NAME(ring-bonds);
CL_DEFUN core::List_sp RingFinder_O::ringBonds(core::List_sp atoms) {
  core::List_sp ringBonds = _Nil<core::T_O>();
  gctools::SmallOrderedSet<Atom_sp> atomSet;
  {
    _BLOCK_TRACE(BF("Put atoms into set"));
    for (auto cur : atoms) {
      Atom_sp atom = oCar(cur).as<Atom_O>();
      LOG(BF("Atom in ring: %s") % atom->description());
      atomSet.insert(atom);
    }
  }
  int num = 0;
  {
    _BLOCK_TRACE(BF("Loop over all bonds and find those that join ring atoms"));
    for (auto cur : atoms) {
      Atom_sp atom = oCar(cur).as<Atom_O>();
      VectorBond bonds = atom->getBonds();
      for (gctools::Vec0<Bond_sp>::const_iterator bi = bonds.begin(); bi != bonds.end(); bi++) {
        LOG(BF("Testing %s") % (*bi)->describeOther(atom));
        if (atomSet.count((*bi)->getAtom1()) && atomSet.count((*bi)->getAtom2())) {
          LOG(BF("It's part of the ring!!!"));
          core::Cons_sp one = core::Cons_O::create((*bi), ringBonds);
          ringBonds = one;
          num++;
        }
      }
    }
  }
  LOG(BF("There were %d bonds in the ring") % num);
  return ringBonds;
}



};
