/*
    File: macroModel.cc
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
#define DEBUG_LEVEL_NONE       
//
// (C) 2004 Christian E. Schafmeister
//

/*
 *	macroModel.cc
 *
 *	Functions that read and write macroModel files.
 *
 */

#include <stdio.h>
#include <string.h>
#include <string>
#include <map>
#include <iostream>
#include <clasp/core/common.h>
#include <cando/chem/matter.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/chem/macroModel.h>
#include <clasp/core/wrappers.h>

namespace chem
{


typedef	struct	{
	int		iAtom;
	int		iBondOrder; //	BondOrder	iBondOrder;
} MMBondType;

typedef	struct	{
	Atom_sp		atom;
	int		iName;
	int		iType;
	MMBondType	bonds[6];
	double		x,y,z;
	int		resNum;
	char		resCode;
	char		resChain;
	int		atomColor;
	double		dCharge;
	double		dMultiPoleCharge;
	char		caResName[16];
	char		caAtomName[16];
} MMAtomType;


typedef	vector<MMAtomType>	MMAtoms;










//
//	openRead
//
//	Open the file for reading
//
void	MacroModelFile_O::openRead( const string& fn )
{
char	caLine[1024];
int	iAtoms;
    this->file = fopen(fn.c_str(),"r");
    fgets( caLine, sizeof(caLine), this->file );
    sscanf( caLine, "%d", &iAtoms );
    this->lines = iAtoms;
}




//
//	openWrite
//
//	Open the file for writing
//
void	MacroModelFile_O::openWrite(const string& fn )
{
    this->file = fopen(fn.c_str(),"w");
}



//
//	atEnd
//
//	Return true if we are at the end of the file
bool	MacroModelFile_O::atEnd()
{
    if ( this->file ) {
        return feof(this->file);
    } else {
	return false;
    }
}


//
//	close
//
void	MacroModelFile_O::close()
{
    fclose(this->file);
    this->file = NULL;
}







//
//	readMolecule( )
//
//
Molecule_sp	MacroModelFile_O::readMolecule()
{
MMAtoms			atoms;
MMAtomType		oneAtom;
Atom_sp			aPAtom;
char			caLine[1024];
string			sTemp;
Vector3			vPos;
int			i,j;
int			iAtoms;
int			iTo;
BondOrder		iOrder;
Atom_sp			aPTo;
gctools::SmallMap<int,Residue_sp>	residues;
map<int,string>		mmodToAmberTypes;
map<int,string>		mmodToAmberElement;
Molecule_sp		mol;
Residue_sp		res;
int			offset, len;
char			caPart[254], *cPTemp;

    mmodToAmberTypes[2] = "C";
    mmodToAmberTypes[3] = "CT";
    mmodToAmberTypes[4] = "CH";
    mmodToAmberTypes[5] = "C2";
    mmodToAmberTypes[6] = "C3";
    mmodToAmberTypes[7] = "CD";
    mmodToAmberTypes[15]= "O";
    mmodToAmberTypes[16]= "OH";
    mmodToAmberTypes[18]= "OC";
    mmodToAmberTypes[25]= "N";
    mmodToAmberTypes[26]= "NT";
    mmodToAmberTypes[31]= "NA";
    mmodToAmberTypes[32]= "N3";
    mmodToAmberTypes[41]= "HC";
    mmodToAmberTypes[42]= "HO";
    mmodToAmberTypes[43]= "H";
    mmodToAmberTypes[44]= "H3";
    mmodToAmberTypes[49]= "S";
    mmodToAmberTypes[53]= "P";
    mmodToAmberTypes[63]= "LP";

    mmodToAmberElement[2] = "C";
    mmodToAmberElement[3] = "C";
    mmodToAmberElement[4] = "C";
    mmodToAmberElement[5] = "C";
    mmodToAmberElement[6] = "C";
    mmodToAmberElement[7] = "C";
    mmodToAmberElement[15]= "O";
    mmodToAmberElement[16]= "O";
    mmodToAmberElement[18]= "O";
    mmodToAmberElement[25]= "N";
    mmodToAmberElement[26]= "N";
    mmodToAmberElement[31]= "N";
    mmodToAmberElement[32]= "N";
    mmodToAmberElement[41]= "H";
    mmodToAmberElement[42]= "H";
    mmodToAmberElement[43]= "H";
    mmodToAmberElement[44]= "H";
    mmodToAmberElement[49]= "S";
    mmodToAmberElement[53]= "P";
    mmodToAmberElement[63]= "LP";


//
// _ENTRY HEADER LINE_ (first line of file or entry)
//
//    FORMAT(1X,I5,2X,A70)
//
//    I5 NIATMS (INTEGER*4) - The number of atoms in the entry
//    A70 (CHARACTER*70) - The text name of the entry (optional)
//
// _ATOM ENTRIES_
// Atom entries start from line 2 and continue until line
//	NIATMS + 1. Each line represents one atom.
//
// FORMAT(1X,I3,6(1X,I5,1X,I1),1X,3(F11.6,1X),I5,2A1,I4,2F9.5,1X,A4,1X,A4)
//
//    I3 (INTEGER*4) The MacroModel atom type number (e.g. 3 = sp3 carbon)
//    6(1X,I5,1X,I1) Describing up to 6 attached atoms as follows:
//    I5 (INTEGER*4) The atom sequence number of an attached (bonded) atom
//    I1 (INTEGER*4) The order of the bond to the attached atom
//		(0 = zero order [used for transition states and
//		coordination complexes], 1 = single bond, 2 = double bond,
//		3 = triple bond).
//    3(F11.6,1X ) Cartesian coordinates for this atom:
//    F11.6 (REAL*4) - X-coordinate in Angstroms.
//    F11.6 (REAL*4) - Y-coordinate in Angstroms.
//    F11.6 (REAL*4) - Z-coordinate in Angstroms.
//    I5 (INTEGER*4) - Optional residue number (used only with biopolymers).
//    A1 (CHARACTER*1) - Optional single-character MacroModel residue
//		code (see below).
//    A1 (CHARACTER*1) - Optional single-character residue chain name.
//    I4 (INTEGER*4) - Optional atom color (see below). If no color
//		is given, then MacroModel will give the atom a standard
//		atom type-based color.
//    2F9.5 (REAL*4) - Optional atomic charge in electron units (BatchMin
//		will recompute charges unless CHGF option used). First
//		charge is used for standard charge/charge electrostatics,
//		second charge is used in conjunction with charge/multipole
//		electrostatics.
//    A4 (CHARACTER*4) - Optional ascii residue name (e.g. PDB residue name)
//    A4 (CHARACTER*4) - Optional ascii atom name (e.g. PDB atom name)
//

    iAtoms = this->lines;
    for ( i=0; i< iAtoms; i++ ) {
	fgets( caLine, sizeof(caLine), this->file );
        LOG("Raw line: %s" , (caLine ) );
//	sscanf( caLine, " %3d %5d %1d %5d %1d %5d %1d %5d %1d %5d %1d %5d %1d %11.6lf %11.6lf %11.6lf %5d%c%c%4d%9.5lf%9.5lf %4s %4s\n",
	sscanf( caLine, " %3d %5d %1d %5d %1d %5d %1d %5d %1d %5d %1d %5d %1d",
			&oneAtom.iType,
			&oneAtom.bonds[0].iAtom, &(oneAtom.bonds[0].iBondOrder),
			&oneAtom.bonds[1].iAtom, &(oneAtom.bonds[1].iBondOrder),
			&oneAtom.bonds[2].iAtom, &(oneAtom.bonds[2].iBondOrder),
			&oneAtom.bonds[3].iAtom, &(oneAtom.bonds[3].iBondOrder),
			&oneAtom.bonds[4].iAtom, &(oneAtom.bonds[4].iBondOrder),
			&oneAtom.bonds[5].iAtom, &(oneAtom.bonds[5].iBondOrder) );

	offset = 53;
	len = 11;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.x = atof(caPart);
	offset += len+1;

	len = 11;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.y = atof(caPart);
	offset += len+1;

	len = 11;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.z = atof(caPart);
	offset += len+1;

	len = 5;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.resNum = atoi(caPart);
	offset += len;

	len = 1;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.resCode = caPart[0];
	offset += len;

	len = 1;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.resChain = caPart[0];
	offset += len;

	len = 4;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.atomColor = atoi(caPart);
	offset += len;

	len = 9;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.dCharge= atof(caPart);
	offset += len;

	len = 9;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	oneAtom.dMultiPoleCharge= atof(caPart);
	offset += len+1;

	len = 4;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	cPTemp = caPart;
	while ( *cPTemp != ' ' && *cPTemp != '\0' ) cPTemp++;
	*cPTemp = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	strcpy( oneAtom.caResName, caPart );
	offset += len+1;

	len = 4;
	strncpy( caPart, &caLine[offset], len );
	caPart[len] = '\0';
	cPTemp = caPart;
	while ( *cPTemp != ' ' && *cPTemp != '\0' ) cPTemp++;
	*cPTemp = '\0';
	LOG("LEFT: |%s|" , (caPart ) );
	strcpy( oneAtom.caAtomName, caPart );
	offset += len+1;

#if 0
			&oneAtom.resNum,	// %5d
			&oneAtom.resCode,	// %c
			&oneAtom.resChain,	// %c
			&oneAtom.atomColor,	// %4d
			&oneAtom.dCharge,	// %9.5lf
			&oneAtom.dMultiPoleCharge,	// %9.5lf
			oneAtom.caResName,	// %4s
			oneAtom.caAtomName	// %4s
#endif
	oneAtom.iName = i+1;
	atoms.push_back(oneAtom);
    }

     LOG("proc.   : %3d %5d %1d %5d %1d %5d %1d %5d %1d %5d %1d %5d %1d %11.6lf %11.6lf %11.6lf %5d%c%c%4d%9.5lf%9.5lf %4s %4s\n"
     , oneAtom.iType %
     oneAtom.bonds[0].iAtom% oneAtom.bonds[0].iBondOrder%
			oneAtom.bonds[1].iAtom% oneAtom.bonds[1].iBondOrder%
			oneAtom.bonds[2].iAtom% oneAtom.bonds[2].iBondOrder%
			oneAtom.bonds[3].iAtom% oneAtom.bonds[3].iBondOrder%
			oneAtom.bonds[4].iAtom% oneAtom.bonds[4].iBondOrder%
			oneAtom.bonds[5].iAtom% oneAtom.bonds[5].iBondOrder%
			oneAtom.x% oneAtom.y% oneAtom.z%
			oneAtom.resNum%	// %5d
			oneAtom.resCode%	// %c
			oneAtom.resChain%	// %c
			oneAtom.atomColor%	// %4d
			oneAtom.dCharge%	// %9.5lf
			oneAtom.dMultiPoleCharge%	// %9.5lf
			oneAtom.caResName%	// %4s
			oneAtom.caAtomName	// %4s
	 );

     mol = Molecule_O::create();
    LOG("Created molecule" );
		//
		// Now create the residue for the atoms
		//
    for ( i=0; i<iAtoms; i++ ) {
	if (residues.count(atoms[i].resNum)==0) {
	    res = Residue_O::create();
	    sTemp = atoms[i].caResName;
	    res->setName(chemkw_intern(sTemp));
	    residues[atoms[i].resNum] = res;
	    mol->addMatter(res);
	}
    }
    LOG("Created residues" );

		//
		// Now create the atoms and stick them in their
		// respective residues
		//

    for ( i=0; i<iAtoms; i++ ) {
	LOG("Creating atom: %d" , (i ) );
	aPAtom = Atom_O::create();
	LOG("Creating atom" );
	if ( residues.count(atoms[i].resNum) > 0 ) {
	    res = residues[atoms[i].resNum];
	} else {
	    LOG("Bad RESIDUE" );
	    exit(1);
	}
	res->addAtom(aPAtom);
	atoms[i].atom = aPAtom;
//	sprintf( caString, "%d", atoms[i].iName );
//	sTemp = caString;
	sTemp = atoms[i].caAtomName;
	aPAtom->setName(chemkw_intern(sTemp));
	aPAtom->setType(chemkw_intern(mmodToAmberTypes[atoms[i].iType]));
	aPAtom->setElementFromString(mmodToAmberElement[atoms[i].iType]);
	vPos = Vector3( atoms[i].x, atoms[i].y, atoms[i].z );
	aPAtom->setPosition(vPos);
        aPAtom->setCharge(atoms[i].dCharge );
    }
    LOG("Created atoms" );

		// Connect bonds
    for ( i=0; i<iAtoms; i++ ) {
	aPAtom = atoms[i].atom;
	for ( j=0; j<6; j++ ) {
	    if ( atoms[i].bonds[j].iAtom != 0 ) {
		iTo = atoms[i].bonds[j].iAtom;
		iOrder = (BondOrder)(atoms[i].bonds[j].iBondOrder);
		if ( iTo >= 1 ) {
		    aPTo = atoms[iTo-1].atom;
		    if ( aPAtom->atLowerUniqueAtomOrderThan(aPTo) )
		    {
			aPAtom->bondTo(aPTo,iOrder);
		    }
		}
	    }
	}
    }
    LOG("Created bonds" );

		// Read the number of lines in the next molecule

    fgets( caLine, sizeof(caLine), this->file );
    sscanf( caLine, "%d", &iAtoms );
    this->lines = iAtoms;
    return mol;
}




};
