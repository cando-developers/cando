       
       
//
// (C) 2004 Christian E. Schafmeister
//




#ifndef	TORSION_DRIVER_H
#define	TORSION_DRIVER_H

#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <exception>
#include "core/common.h"
#include "core/numerics.h"
//#include "matter.h"
//#include "coordSys.h"

#include "chemPackage.h"


namespace chem {

SMART(Matter);
SMART(Atom);
SMART(CoordinateSystem);



class	TorsionAtom {
public:
	Atom_sp	getAtom() { return this->atom; };
	void	setAtom(Atom_sp a) { this->atom = a; };
	Vector3	getUntransformed() { return this->untransformed; };
	void	setUntransformed(double x, double y, double z) {
		this->untransformed.set(x,y,z);};
	void	setUntransformed(const Vector3& a) {
		this->untransformed=a;};

private:
	Atom_sp			atom;
	Vector3			untransformed;
};

typedef	vector<TorsionAtom>::iterator	TorsionAtomIterator;

class	TorsionSegment {
public:
	Atom_sp				fixed;
	Atom_sp				movable;
	BondOrder			order;
	int				steps;
	vector<TorsionAtom>		atoms;	//!< Atoms to rotate
	CoordinateSystem_sp		coord;
	Matrix				relative;
	vector<Matrix>			rotations; //!< A Matrix for every angle that we rotate this segment by
	vector<Matrix>::iterator	rotation;
	Matrix				accumulated;
			//
			// Define a matrix update tree
			//
	int				parent;
	vector<int>			children;
};






class	Dumb_TorsionDriver;
    typedef	gctools::smart_ptr<Dumb_TorsionDriver>	RPTorsionDriver;

class	Dumb_TorsionDriver{
private:
	bool			driveMode;
	Aggregate_sp		aggregate;
	vector<TorsionSegment>	torsions;

	void insertTorsion(	Atom_sp	fixed,
				Atom_sp	movable,
				int	steps );
    	core::LongLongInt basicCountConformations(int segment);
	int	countOfAtomsWithName( Aggregate_sp agg, const string& name );
	Atom_sp	firstAtomWithName( Aggregate_sp agg, const string& name );
	void	updateAccumulatedTransforms(int start);
	void	calculateNewCoordinates(int start);

public:

    friend	RPTorsionDriver	TorsionDriver();

    void	setAggregateAndRoot( Aggregate_sp	agg,
					string	rootName );

    void	addTorsion(	string	fixed,
				string	movable,
				int	steps );

    int		numberOfSegments() { return this->torsions.size(); };
    int		firstValidSegment() { return 1; };
    Atom_sp	fixedAtomOfSegment(int s) { return this->torsions[s].fixed;};
    Atom_sp	movableAtomOfSegment(int s) { return this->torsions[s].movable;};

    int		numberOfAtomsInSegment( int segment ) { return this->torsions[segment].atoms.size();};
    Atom_sp	atomInSegment(int segment, int iAtom ) { return this->torsions[segment].atoms[iAtom].getAtom();};

    void	dump();
    string	dumpString();

    void	prepareToDriveTorsions();

//    addon::QDomNode_sp	renderXml();


    core::LongLongInt countConformations() { return this->basicCountConformations(1);};


		//
		// Advance to the next conformation
		//
		// This uses an odometer algorithm and returns
		// the lowest segment index that was updated.
		//
		// If tests show that a particular segment is colliding
		// with something that it shouldn't then call
		// advanceToNextConformation with the offending segment
		// and it will be incremented.  This is a sort of
		// dead-end-elimination.
		//
		//
    int		highestSegmentIndex() { return this->torsions.size()-1; };
    int		advanceToNextConformationStartFromSegment(int segment);
    int		advanceToNextConformation() { return this->advanceToNextConformationStartFromSegment(this->highestSegmentIndex());};

    TorsionAtomIterator	atomBegin(int segment) { return this->torsions[segment].atoms.begin(); };
    TorsionAtomIterator	atomEnd(int segment) { return this->torsions[segment].atoms.end(); };

    Dumb_TorsionDriver() {
		REF_CREATE("Dumb_TorsionDriver");
		this->driveMode = false;
    };
    virtual ~Dumb_TorsionDriver(){REF_DELETE("Dumb_TorsionDriver");};

};




inline	RPTorsionDriver	TorsionDriver() {
    return RPTorsionDriver(new Dumb_TorsionDriver());
}

};
#endif

