/*
    File: torsionDriver.cc
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


#define	TURN_DEBUG_OFF

#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <exception>
#include <core/boostSmartPointers.h>
#include <cando/chem/matter.h>
#include <clasp/core/exceptions.h>
#include <cando/chem/loop.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/torsionDriver.h>
#include <cando/chem/coordSys.h>

using namespace chem;


int	Dumb_TorsionDriver::countOfAtomsWithName( Aggregate_sp	agg,
							const string& name )
{
Loop	lAtoms;
Atom_sp	a;
int	count;

    count = 0;
    lAtoms.loopTopGoal(agg,ATOMS);
    while ( lAtoms.advanceLoopAndProcess() ) {
	a = lAtoms.getAtom();
	if ( a->getName() == name ) count++;
    }
    return count;
}


#if 0 // DISCARD
Atom_sp	Dumb_TorsionDriver::firstAtomWithName( Aggregate_sp	agg,
						const string& name )
{
Loop	lAtoms;
Atom_sp	a;
int	count;
bool	gotAtom;
    gotAtom = false;
    count = 0;
    lAtoms.loopTopGoal(agg,ATOMS);
    while ( lAtoms.advanceLoopAndProcess() ) {
	a = lAtoms.getAtom();
	if ( a->getName() == name ) {
	    gotAtom = true;
	    break;
	}
    }
    if ( !gotAtom ) {
	SIMPLE_ERROR(("Could not find atom with name: ", name ));
    }
    return a;
}
#endif //DISCARD









//
//	insertTorsion
//
//	Insert the torsion into the TorsionSegment vector.
//
//	The TorsionSegments are inserted in order of how far they
//	fall from the root.  This is indicated by the value of each atoms
//	(tempInt) variable. This was set by a spanningLoop that started
//	on the root atom. By inserting the TorsionSegment in ascending order
//	by the (tempInt) variable, they will be in the correct order for
//	torsion driving.
//
void	Dumb_TorsionDriver::insertTorsion(	Atom_sp	fixed,
						Atom_sp	movable,
						int	steps )
{
vector<TorsionSegment>::iterator	pos;
TorsionSegment				torsion;

	//
	// Skip the first torsion, it's a dummy
	//
    for ( pos=this->torsions.begin()+1;
		pos!=this->torsions.end();
		pos++ ) {
	if ( fixed->getTempInt() < pos->fixed->getTempInt() ) break;
    }
    torsion.fixed = fixed;
    torsion.movable = movable;
    torsion.steps = steps;
    torsion.parent = -1;
    if ( pos==this->torsions.end() ) {
	this->torsions.push_back(torsion);
    } else {
	this->torsions.insert(pos,torsion);
    }
}






void	Dumb_TorsionDriver::updateAccumulatedTransforms(int start)
{
vector<TorsionSegment>::iterator	segment;
vector<int>::iterator	pos;
TorsionSegment		*parent;
Matrix			mat;

    for ( segment=this->torsions.begin()+start;
		segment != this->torsions.end();
		segment++ ) {
	parent = &(this->torsions[segment->parent]);
	mat = segment->relative*(*(segment->rotation));
	segment->accumulated = parent->accumulated*mat;
    }
	//
	// First update my accumulated transform
	//
#if 0 //[
    segment = &(this->torsions[start]);
    parent = &(this->torsions[segment->parent]);
    mat = segment->relative*(*(segment->rotation));
    segment->accumulated = parent->accumulated*mat;

    LOG("Updated accumulated transform for segment: {}" , (start ) );
#ifdef	DEBUG_ON //[
    LOG("    rotation = " );
	segment->rotation->dump();
    LOG("    accumulated = " );
	segment->accumulated.dump();
#endif //]
	//
	// Then update the children
	//

    for ( pos=this->torsions[start].children.begin();
	    pos!=this->torsions[start].children.end();
	    pos++ ) {
	this->updateAccumulatedTransforms(*pos);
    }
#endif //]
}




//
//	setAggregateAndRoot
//
//	Set the aggregate and the root of the
void	Dumb_TorsionDriver::setAggregateAndRoot( Aggregate_sp	agg,
						string	rootName )
{
SpanningLoop_sp	spanning;
Atom_sp		root, atom;
Loop		lAtoms;
bool		gotRoot;
TorsionSegment	torsion;
int		spanningIndex;

    this->driveMode = false;	// We are about to start adding torsions
    this->torsions.erase(this->torsions.begin(),this->torsions.end());
    this->aggregate = agg; // Should I use a copy?
    gotRoot = false;
    lAtoms.loopTopGoal( this->aggregate, ATOMS );
    while ( lAtoms.advanceLoopAndProcess() ) {
	atom = lAtoms.getAtom();
	if ( atom->getName() == rootName ) {
	    root = atom;
	    gotRoot = true;
	}
    }
    if ( !gotRoot ) {
	SIMPLE_ERROR("Could not find atom name: {} in aggregate" , rootName );
    }

		//
		// Loop over the atoms in the aggregate using a
		// SpanningLoop.  Set the (tempInt) variable of each
		// atom with a counter that will allow the insertTorsion
		// function to properly order torsions in the TorsionSegment
		// list.
		//

    spanning = agg->lisp()->create<SpanningLoop_O>();
    spanning->setTop(root);
    spanningIndex = 0;
    while ( spanning->advanceLoopAndProcess() ) {
	atom = spanning->getAtom();
	atom->setTempInt(spanningIndex);
	spanningIndex++;
    }

		//
		// Add a dummy torsion segment for seamless updating of
		// matrices
		//
    torsion.steps = 1;
    this->torsions.push_back(torsion);


}




void	Dumb_TorsionDriver::addTorsion(	string	atom1,
					string	atom2,
					int	steps )
{
Loop	lAtoms;
Atom_sp	atom, a1, a2;
Atom_sp	fixed, movable, backSpan;

    if ( this->driveMode ) {
	SIMPLE_ERROR("You can't add torsions once you've started driving torsions");
    }
    ANN(this->aggregate);
    if ( this->aggregate.nilp() )
    {
	SIMPLE_ERROR("You must have defined the aggregate and root before adding torsions to drive.");
    }
    if ( this->countOfAtomsWithName(this->aggregate,atom1 ) != 1 ) {
	SIMPLE_ERROR("There must be a unique atom with the name: {}" , atom1 );
    }
    if ( this->countOfAtomsWithName(this->aggregate,atom2 ) != 1 ) {
	SIMPLE_ERROR("There must be a unique atom with the name: " , atom2 );
    }
    a1 = this->aggregate->firstAtomWithName(atom1);
    a2 = this->aggregate->firstAtomWithName(atom2);
    backSpan = a2->getBackSpan();
    if ( a1 == a2->getBackSpan() ) {
	this->insertTorsion( a1, a2, steps );
	return;
    }
    if ( a2 == a1->getBackSpan() ) {
	this->insertTorsion( a2, a1, steps );
    }
    SIMPLE_ERROR("Movable atom has to point back to the fixed atom");
}




//
//	prepareToDriveTorsions
//
//	Prepare this TorsionDriver to drive the torsions that
//	have been specified.
//
void	Dumb_TorsionDriver::prepareToDriveTorsions()
{
vector<TorsionSegment>::iterator	pos, parent;
Vector3					origin;
Vector3					vTemp, xDir, yDir, zDir;
BondOrder				bo;
SpanningLoop_sp				spanning;
Atom_sp					a;
int					index, i;
Matrix					rot;
double					angle;
TorsionAtom				tatom;
vector<TorsionAtom>::iterator		tai;
Matrix					toCanonical;
CoordinateSystem_sp			coord;

    if ( this->driveMode ) {
	SIMPLE_ERROR("You can't prepare to drive torsions when you're already driving them");
    }
	//
	// Break the bonds around which we are rotating the torsions
	//
	// Skip the first segment, it's a dummy
    LOG("Breaking torsion bonds" );
    for ( pos=this->torsions.begin()+1;
		pos!=this->torsions.end();
		pos++ ) {
	bo = pos->movable->bondOrderTo(pos->fixed);
	if ( bo==noBond ) {
	    SIMPLE_ERROR("The two atoms of a torsion must be bonded. Atom: {} is not bonded to {}" , pos->fixed->getName() , pos->movable->getName() );
	}
	pos->order = bo;
	pos->fixed->removeBondBetween(pos->movable,pos->fixed);
    }

	//
	// Now find all the atoms that belong to each segment.
	//
	// Do this by building a SpanningLoop from each movable atom.
	// Now that all the bonds are removed across torsions, each
	// spanning loop will only access the atoms of each segment
	//
	// Set the tempInt variable of each atom to the index of the
	// segment that the atom belongs to. We will use that to
	// define the matrix update tree.
	//
	// Skip the first segment, it's a dummy

    LOG("Identifying atoms that belong to each segment" );
    spanning = SpanningLoop_O::create();
    for ( pos=this->torsions.begin()+1;
		pos!=this->torsions.end();
		pos++ ) {
	spanning->setTop(pos->movable);
	while ( spanning->advanceLoopAndProcess() ) {
	    a = spanning->getAtom();
	    tatom.setAtom(a);
	    tatom.setUntransformed(0,0,0);
	    pos->atoms.push_back(tatom);
	    a->setTempInt(pos-this->torsions.begin());
	}
    }

	//
	// Construct the matrix update tree.
	// First define the parent of each segment.
	//
    LOG("Constructing the matrix update tree - parents" );
    for ( pos=this->torsions.begin(), index = 0;
		pos!=this->torsions.end();
		pos++, index++ ) {
		// Define the index of the parent of this segment
		// That will be -1 if the segment is the first one
		// otherwise its the segment of the fixed atom
	if ( pos == this->torsions.begin() ) {
	    pos->parent = -1;
	} else if ( pos== this->torsions.begin()+1 ) {
	    pos->parent = 0;
	} else {
	    pos->parent = pos->fixed->getTempInt();
	}
    }
	//
	// Now define the children
	//
    LOG("Constructing the matrix update tree - children" );
    for ( pos=this->torsions.begin()+1, index = 1;
		pos!=this->torsions.end();
		pos++, index++ ) {
		// Define the index of the parent of this segment
		// That will be -1 if the segment is the first one
		// otherwise its the segment of the fixed atom
	this->torsions[pos->parent].children.push_back(index);
    }

	//
	// Build a coordinate system for every fixed atom
	// with its Z-axis along the line joining the fixed atom to the
	// movable one and its X-axis perpendicular to the Z-axis.
	//
    LOG("Building a coordinate system for every fixed atom" );
    for ( pos=this->torsions.begin();
		pos!=this->torsions.end();
		pos++ ) {
		//
		// Calculate the zDir and xDir vectors for
		// a coordinate system.
		//
	if ( pos==this->torsions.begin() ) {
	    origin.set(0,0,0);
	    xDir.set(1,0,0);
	    zDir.set(0,0,1);
	    pos->coord = CoordinateSystem_O::create();
	    pos->coord->defineForVectorsOriginXDirZDir(origin,xDir,zDir);
	} else {

	    origin = pos->fixed->getPosition();
	    zDir = pos->movable->getPosition()-pos->fixed->getPosition();
	    zDir = zDir.normalized(_lisp);
	    if ( pos->movable->numberOfBonds() == 0 ) {
		SIMPLE_ERROR("The torsion around {} and {} won't roate anything. Please remove it" 
					% pos->fixed->getName() % pos->movable->getName() );
	    }
	    vTemp = pos->movable->bondedNeighbor(0)->getPosition()
			    - pos->movable->getPosition();
	    xDir = vTemp.crossProduct(zDir);
	    xDir = xDir.normalized(_lisp);

		    //
		    // Create a coordinate system
		    //
	    pos->coord = CoordinateSystem_O::create();
	    pos->coord->defineForVectorsOriginXDirZDir( origin, xDir, zDir );

		    //
		    // Now reconnect the bond
		    //
	    pos->fixed->bondTo( pos->movable, pos->order );
	}
    }

	//
	// Now calculate the relative transformation
	//
    LOG("Calculating the relative transformations" );
    for ( pos=this->torsions.begin()+1;
		pos!=this->torsions.end();
		pos++ ) {
	LOG("  for pos={}" , (pos-this->torsions.begin() ) );
	LOG("  Parent ={}" , (pos->parent ) );
	parent = this->torsions.begin()+pos->parent;
	LOG("  Calculating transform" );
	coord = parent->coord;
	pos->relative = coord->matrixForTransformTo(pos->coord);
	LOG("  Calculating transform done" );
    }

	//
	// Generate the rotation matrices that we will need to iterate over
	//

    LOG("Generating rotation matrices" );
    for ( pos=this->torsions.begin();
		pos!=this->torsions.end();
		pos++ ) {
	LOG("Generating rotation matrices for segment: {}" , (pos-this->torsions.begin() ) );
	for ( i=0; i<pos->steps; i++ ) {
	    angle = i*(0.0174533*360.0/pos->steps);
	    LOG("Rotating by {} radians" , (angle ) );
	    rot.rightHandedRotationZ(-angle);
	    pos->rotations.push_back(rot);
	}
        pos->rotation = pos->rotations.begin();
    }

	//
	// Update the accumulated transforms.
	// The first accumulated matrix gets us from the origin
	// to the first coordinate system.

    LOG("Updating the accumulated transformations" );
		// Get the first real segment
    pos = this->torsions.begin()+1;
    this->torsions.begin()->accumulated.setToIdentity(); // = pos->coord->matrixFromCanonical();

    LOG("Updating the accumulated transformations" );
    this->updateAccumulatedTransforms(1);

	//
	// Now calculate positions of segment atoms transformed into
	// the cannonical coordinate system
	//

    LOG("Calculating the positions of segment atoms in the cannonical coord system" );
    for ( pos=this->torsions.begin()+1; pos!=this->torsions.end(); pos++ ) {
	toCanonical = pos->coord->matrixToCanonical();
	for ( tai=pos->atoms.begin(); tai!=pos->atoms.end(); tai++ ) {
	    tai->setUntransformed(toCanonical*(tai->getAtom()->getPosition()));
	}
    }
    LOG("Ready to drive torsions" );
    this->driveMode = true;
}





#if 0 //[
//
//	renderXml
//
//	Render a representation of the torsions that we are going to
//	drive for debugging.
//
//
adapt::QDomNode_sp	Dumb_TorsionDriver::renderXml()
{
adapt::QDomNode_sp	graphics, aggGraphics;
adapt::QDomNode_sp	coords;
adapt::QDomNode_sp	gr;
vector<TorsionSegment>::iterator	pos;
Vector3		origin, pnt;
vector<TorsionAtom>::iterator		tai;


    if ( !this->driveMode ) {
	SIMPLE_ERROR("TorsionDriver must be in drive mode to render graphics"));
    }

    graphics = xmlGraphics("torsions");

		//
		// Render the aggregate
		//

    aggGraphics = this->aggregate->asXmlWithCoordinates();
    graphics->addChild(aggGraphics);


		//
		// For each torsion render a coordinate system
		//

    for ( pos=this->torsions.begin(); pos!=this->torsions.end(); pos++ ) {
	gr = pos->coord->renderXml(1.0,1.0);
	graphics->addChild(gr);
	for ( tai=pos->atoms.begin(); tai!=pos->atoms.end(); tai++ ) {
	    gr = xmlSphere(XML_YELLOW,tai->getUntransformed(),0.3);
	    graphics->addChild(gr);
	}
    }

    coords = xmlBlock();
    origin.set(0,0,0);
    pnt.set(10,0,0);
    gr = xmlLine( XML_RED, 2, origin, pnt );
    coords->addChild(gr);
    pnt.set(0,20,0);
    gr = xmlLine( XML_GREEN, 2, origin, pnt );
    coords->addChild(gr);
    pnt.set(0,0,30);
    gr = xmlLine( XML_BLUE, 2, origin, pnt );
    coords->addChild(gr);
    graphics->addChild(coords);
    return graphics;
}
#endif //]


//
//	calculateNewCoordinates
//
//	Update the coordinates for the atoms
//
void	Dumb_TorsionDriver::calculateNewCoordinates(int start)
{
vector<TorsionSegment>::iterator	pos;
vector<TorsionAtom>::iterator		tai;
Vector3					trans;

    for ( pos=this->torsions.begin()+start; pos!=this->torsions.end(); pos++ ) {
	LOG("Applying matrix at segment: {}" , (pos-this->torsions.begin() ) );
#ifdef	DEBUG_ON
	pos->accumulated.dump();
#endif
	for ( tai=pos->atoms.begin(); tai!=pos->atoms.end(); tai++ ) {
	    LOG("Transforming atom: {}" , (tai->getAtom()->getName().c_str() ) );
	    trans = pos->accumulated*(tai->getUntransformed());
	    tai->getAtom()->setPosition(trans);
	}
    }
}


//
//	advanceToNextConformationStartFromSegment
//
//	Advance to the next conformation.
//	Start at the specified segment and work backwards.
//	Return the index of the most significant segment to be updated.
//	A return value of 0 means that no more conformations will be generated.
//
int	Dumb_TorsionDriver::advanceToNextConformationStartFromSegment(int segment)
{
vector<TorsionSegment>::iterator	pos;
vector<int>::iterator			cpos;

    LOG("----Advancing from segment: {} highest: {}" , (segment) , (this->torsions.end()-this->torsions.begin()-1) );
    if (segment > this->highestSegmentIndex() ) {
	segment = this->highestSegmentIndex();
	LOG("Trimmed to start at segment: {}" , (segment ) );
    } else if ( segment < this->highestSegmentIndex()) {
	LOG("Resetting rotations for segments after {}" , (segment ) );
	for ( pos=this->torsions.begin()+segment+1;
		pos!=this->torsions.end();
		pos++ ) {
	    pos->rotation = pos->rotations.begin();
	}
#if 0		// zero children only
        for ( cpos=this->torsions[segment].children.begin();
		cpos!=this->torsions[segment].children.end();
		cpos++ ) {
	    this->resetRotation(*cpos);
	}
#endif
    }

    pos = this->torsions.begin()+segment; // this->torsions.end()-1;
    while ( pos > this->torsions.begin() ) {
	LOG("Incrementing segment: {}" , (pos-this->torsions.begin() ) );
	// Increment the current torsionSegment rotation
	pos->rotation++;
	LOG("Incremented rotation {} of {}" , (pos->rotation-pos->rotations.begin()) , (pos->rotations.end()-pos->rotations.begin() ) );
	if ( pos->rotation != pos->rotations.end() ) break;
	LOG("Carried" );
	pos->rotation = pos->rotations.begin();
	pos--;
    }
    if ( pos > this->torsions.begin() ) {
	LOG("Updating transforms and coordinates from segment: {}" , (pos-this->torsions.begin() ) );
	this->updateAccumulatedTransforms(pos-this->torsions.begin());
	this->calculateNewCoordinates(pos-this->torsions.begin());
	return pos-this->torsions.begin();
    }
    return 0;
}




//
//	dump
//
//	Dump a representation of this torsion driver
//
void	Dumb_TorsionDriver::dump()
{
vector<TorsionSegment>::iterator	pos;
vector<TorsionAtom>::iterator		ap;
vector<int>::iterator			ip;
vector<Matrix>::iterator		mp;

    if ( this->driveMode ) {
      core::clasp_write_string("Torsion driver in DRIVE mode\n");
    } else {
      core::clasp_write_string("Torsion driver in ADD_TORSION mode\n");
    }
    for ( pos=this->torsions.begin(); pos != this->torsions.end(); pos++ ) {
//	printf( "Driven torsion #%lX\n", (pos-this->torsions.begin() );
	if ( pos!=this->torsions.begin() ) {
          core::clasp_write_string(fmt::format( "  Parent: {}\n" , pos->parent ));
          core::clasp_write_string(fmt::format( "      Fixed: {}\n" , pos->fixed->getName().c_str() ));
          core::clasp_write_string(fmt::format( "    Movable: {}\n" , pos->movable->getName().c_str() ));
          core::clasp_write_string(fmt::format( "      Steps: {}\n" , pos->steps ));
          core::clasp_write_string("      Atoms:\n");
          for ( ap=pos->atoms.begin(); ap!=pos->atoms.end(); ap++ ) {
            core::clasp_write_string(fmt::format(" {}\n", ap->getAtom()->getName() ));
          }
          core::cl__terpri();
          core::clasp_write_string("      Children:\n");
          for ( ip=pos->children.begin(); ip!=pos->children.end(); ip++ ) {
            core::clasp_write_string(fmt::format(" {:3}\n", *ip ));
          }
          core::cl__terpri();
          core::clasp_write_string("      coordinate system:\n");
          pos->coord->dump();
          core::clasp_write_string("      relative matrix:\n");
          pos->relative.dump();
          for (mp=pos->rotations.begin(); mp!=pos->rotations.end(); mp++ ) {
            mp->dump();
          }
	}
        core::clasp_write_string("      accumulated matrix:\n");
	pos->accumulated.dump();
        core::cl__terpri();
    }
}

//
//	dumpString
//
//	Dump a representation of this torsion driver
//
string	Dumb_TorsionDriver::dumpString()
{
vector<TorsionSegment>::iterator	pos;
vector<TorsionAtom>::iterator		ap;
vector<int>::iterator			ip;
stringstream				ss;

    ss.str("rots: ");
    for ( pos=this->torsions.begin()+1; pos < this->torsions.end(); pos++ ) {
	ss << pos->fixed->getName();
	ss << "-";
	ss << pos->movable->getName();
	ss << "(";
	ss << (pos->rotation-pos->rotations.begin()) << ") ";

    }
    ss << "\n";
    return ss.str();
}





//
//	basicCountConformations
//
//	Return the number of conformations that will be generated
//
core::LongLongInt Dumb_TorsionDriver::basicCountConformations(int segment)
{
core::LongLongInt 		count;
vector<int>::iterator	child;

    count = this->torsions[segment].steps;
    for ( child=this->torsions[segment].children.begin();
		child!=this->torsions[segment].children.end();
		child++ ) {
	count *= this->basicCountConformations(*child);
    }
    return count;
}
