/*
    File: octree.h

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
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */


# ifndef OCTREE_H
# define OCTREE_H


#include <clasp/core/common.h>
#include <cando/chem/atom.h>
#include <cando/geom/vector3.h>
#include <cando/chem/chemPackage.h>


namespace chem {

  FORWARD(OctNode);
  class OctNode_O : public core::CxxObject_O {
    LISP_CLASS(chem,ChemPkg,OctNode_O,"OctNode",core::CxxObject_O);
  public:
    bool fieldsp() const { return true;};
    void fields(core::Record_sp record);
  public:
    int                         iNodeNum;
    int                         iStatus;
    int                         iDepth;
    Vector3                     vCorner;
    int                         iAtoms;
    gctools::Vec0<Atom_sp>      PaAtomList;
    gctools::Vec0<float>        _PfCharges;       /* float to save space */
    OctNode_sp                  PonChildren[8];
  OctNode_O() : iDepth(0), iAtoms(0),
      PonChildren{_Unbound<OctNode_O>(),_Unbound<OctNode_O>(),
        _Unbound<OctNode_O>(),_Unbound<OctNode_O>(),
        _Unbound<OctNode_O>(),_Unbound<OctNode_O>(),
        _Unbound<OctNode_O>(),_Unbound<OctNode_O>()} {};
    string __repr__() const;
  };


  FORWARD(Octree);
  class Octree_O : public core::CxxObject_O {
    LISP_CLASS(chem,ChemPkg,Octree_O,"Octree",core::CxxObject_O);
  public:
    bool fieldsp() const { return true;};
    void fields(core::Record_sp record);
  public:
  int                           iType;  
  int                           _iMaxDepth;
  int                           iTreePoints;
  int                           iDielectric;
  size_t                        iNodeCount;
  float                         fVolume;
  double                        dGridSize;
  float                 	fNewCharge;
  float                         fMinCharge;
  float                         fMaxCharge;
  Vector3                       vMinCharge;
  Vector3                       vMaxCharge;
  Vector3                       vNewPoint;
  Atom_sp                	aClosestAtom;
  double                	dClosestDistance;
  gctools::Vec0<int>            PiDensities;
  gctools::Vec0<int>            iNodeNumCharges;
  gctools::Vec0<float>          _PfCharges;       /* float to save space */
  gctools::Vec0<double>         PdHalfEdges;
  gctools::Vec0<double>         PdHalfDiagonals;
  gctools::Vec0<Atom_sp>        vaAtoms;
  OctNode_sp                    onHead;

  double dDistanceSq( Vector3 Pv1, Vector3 Pv2 );
  double dDistance( Vector3 Pv1, Vector3 Pv2 );
  //static void dumpNode( OctNode_sp PonNode );
  void dumpNode( OctNode_sp PonNode);
  void dumpOctree();
  //static OctNode_sp PonMakeChildren( Vector3& vCorner, int iDepth, int iStatus );
  //static int dFinalCheck( OctNode_sp PonNode, int iAtoms, Atom_sp PaAtomList );
  //static void DestroyOctant( OctNode_sp PonNode, int iStatus );
  //static int iBuildShellOctant( OctNode_sp PonNode, int iAtoms, Atom_sp *PaAtomList );
  //static int iBuildInteriorOctant( OctNode_sp PonNode, int iAtoms, Atom_sp *PaAtomList );
  //OctNode_sp PonMakeChildren( Vector3& vCorner, int iDepth, int iStatus );
  void PonMakeChildren(OctNode_sp ponNode, int iDepth, int iStatus, size_t& iNodeCount );
  int dFinalCheck( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp> PaAtomList, double dShellRadius  );
  void DestroyOctant( OctNode_sp PonNode, int iStatus );
  int iBuildShellOctant( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp> PaAtomList,  double dShellRadius );
  int iBuildInteriorOctant( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp> PaAtomList );
  void OctTreeDestroy();// Octree_sp PoctTree );
  //static void OctNodeInitCharges( OctNode_sp *PonNode );
  //static void OctNodePrintGrid( OctNode_sp *PonNode );
  //static void SplitIncludedNode( OctNode_sp PonNode );
  //static int OctNodeDeleteSphere( OctNode_sp PonNode );
  //static void OctNodeUpdateCharge( OctNode_sp *PonNode, int iParentAtoms, Atom_sp *PaParentAtoms );
  //static int OctNodeCheckSolvent( OctNode_sp PonNode );
  void OctNodeInitCharges( OctNode_sp PonNode, int iDistanceCharge);
  void OctNodePrintGrid( OctNode_sp PonNode, int iColor );
  void SplitIncludedNode( OctNode_sp PonNode );
  int OctNodeDeleteSphere( OctNode_sp PonNode, double dDeleteRadius  );
  void OctNodeUpdateCharge( OctNode_sp PonNode, int iParentAtoms, gctools::Vec0<Atom_sp> PaParentAtoms, int iDistanceCharge );
  int OctNodeCheckSolvent( OctNode_sp PonNode );
  void octOctTreeCreate(Aggregate_sp uUnit, int iType, double dGridSpace, 
			double dAddExtent, double dShellExtent, int iIncludeSolvent, bool bVerbose);
  core::T_mv OctTreeInitCharges(/*Octree_sp octTree,*/ int iAtomOption, int iDielectric, 
                           double dCutDist);	
  void OctTreeDescribe();	/* ( Octree_sp ) */
  void OctTreeDeleteSphere(/*Octree_sp octTree, */ Vector3 vPoint, double dRadius);
  core::T_mv OctTreeUpdateCharge(/*Octree_sp octTree,*/ Vector3 vNewPoint, 
                            float fCharge, double dCutDist);
  Residue_sp rOctTreeCheckSolvent(/*Octree_sp octTree,*/ Vector3 vPoint);
  
//  void OctTreePrintGrid(Octree_sp octTree, core::T_sp stream, int iColor);

  };

/* OCTREEt.iOctType types */

#define OCT_SHELL 		1
#define OCT_INTERIOR_SOLUTE 	2
#define OCT_INTERIOR_SOLVENT	3

/* atom selection options */

#define AT_OCTREE	1

/* dielectric options */

#define DIEL_R2		1
#define DIEL_R		2

/* types of printing */

#define COLOR_CUT	1
#define COLOR_RANGE	2
#define COLOR_DEPTH	3
#define COLOR_NONE	4



}; // namespace chem

TRANSLATE(chem::OctNode_O);
TRANSLATE(chem::Octree_O);

# endif
