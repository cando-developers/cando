

/*
 *      File:   octree.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *
 *  octree.c by Bill Ross
 *
 *	Notes
 *
 *	What is a an octree? Take a cube, divide it into 8 cubes, divide
 *	these into 8, down to the desired limit of resolution. Octrees 
 *	are a simple way to deal with space using minimal resources to 
 *	achieve the desired resolution - the subdivision of a given cube 
 *	can be avoided if the cube can be shown to be homogenous.
 *
 *	Recursive programming can be used to elegantly traverse tree-type 
 *	data structures such as octrees. The basic concept is that a
 *	routine for processing a node calls itself for processing any
 *	sub-nodes. So that recursion does not extend infinitely, a
 *	recursive routine will test some condition or conditions to
 *	determine if the end of recursion is reached. This technique
 *	hasn't been used in Fortran very much because Fortran did not
 *	originally allow multiple calls to a single routine in the
 *	stack, because local variables were not stack-based.
 */

#include <stdio.h>
#include <clasp/core/foundation.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/loop.h>
#include <cando/chem/octree.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/bformat.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/aggregate.h>



namespace chem {

time_t	time_start;

/*	status of octree node	*/
};

#define OCT_UNKNOWN	-1
#define OCT_INCLUDED	1
#define OCT_EXCLUDED	2
#define OCT_PARTIAL	3

SYMBOL_EXPORT_SC_(KeywordPkg,shell);
SYMBOL_EXPORT_SC_(KeywordPkg,interior_solute);
SYMBOL_EXPORT_SC_(KeywordPkg,interior_solvent);
SYMBOL_EXPORT_SC_(ChemPkg,STARoctree_typeSTAR);

CL_BEGIN_ENUM(chem::Octree_O::OctreeType, chem::_sym_STARoctree_typeSTAR, "octree_type");
CL_VALUE_ENUM(kw::_sym_shell,chem::Octree_O::Shell);
CL_VALUE_ENUM(kw::_sym_interior_solute,chem::Octree_O::InteriorSolute);
CL_VALUE_ENUM(kw::_sym_interior_solvent,chem::Octree_O::InteriorSolvent);
CL_END_ENUM(chem::_sym_STARoctree_typeSTAR);
// Define a translator for the enum
#define MY_CL_ENUM_TRANSLATOR(_sym_,_type_) \
namespace translate { \
template <> struct to_object<_type_> \
{								 \
  typedef	_type_	GivenType;	 \
  static core::T_sp convert(const GivenType& val) \
  {_G(); \
    core::SymbolToEnumConverter_sp converter = _sym_->symbolValue().as<core::SymbolToEnumConverter_O>(); \
    return (converter->symbolForEnum(val)); \
  } \
}; \
template <> \
struct from_object<_type_> \
{								 \
  typedef	_type_ 	ExpectedType; \
  typedef	ExpectedType 	DeclareType; \
  DeclareType _v; \
  from_object(gctools::smart_ptr<core::T_O> o) { \
    printf("%s:%d  from_object o -> %s\n", __FILE__, __LINE__, _rep_(o).c_str() ); \
    core::SymbolToEnumConverter_sp converter = _sym_->symbolValue().as<core::SymbolToEnumConverter_O>(); \
    _v = converter->enumForSymbol<_type_>(o.as<core::Symbol_O>()); \
  } \
}; \
};
MY_CL_ENUM_TRANSLATOR(chem::_sym_STARoctree_typeSTAR,chem::Octree_O::OctreeType);


namespace chem{

/* 
 *  globals used for octree evaluation
 */

//double          	dGridSize;
//double          	dShellRadius;
//int             	iColorMethod;

//int             	iTreeGridPoints,
//int                     iNodeCount;
//float           	fVolume;
//gctools::Vec0<int>	PiDensities;
//gctools::Vec0<double> 	PdHalfEdges, PdHalfDiagonals;

//int             	iDistanceCharge;
//int             	iChargeAtoms;
//double          	dCutRadius;
//gctools::Vec0<float>	PfCharges;
//gctools::Vec0<Atom_sp>	PaChargeAtoms;
//float           	fMinCharge, fMaxCharge;
//Vector3	                vMinCharge, vMaxCharge;

//Vector3          	vNewPoint;
//double           	dDeleteRadius;
//float           	fNewCharge;

//Atom_sp          	aClosestAtom;
//double          	dClosestDistance;

//FILE	                *fChargeFile;
//core::T_sp              stream;



void OctNode_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,nodenum), this->iNodeNum);
  node->field( INTERN_(kw,status), this->iStatus);
  node->field( INTERN_(kw,depth), this->iDepth);
  node->field_if_not_default( INTERN_(kw,corner), this->vCorner, Vector3());
  node->field_if_not_default( INTERN_(kw,atoms), this->iAtoms, 0);
  node->field_if_not_empty( INTERN_(kw,atomList), this->PaAtomList);
  node->field_if_not_empty( INTERN_(kw,charges), this->_PfCharges);
  node->field_if_not_unbound(INTERN_(kw,child0),this->PonChildren[0]);
  node->field_if_not_unbound(INTERN_(kw,child1),this->PonChildren[1]);
  node->field_if_not_unbound(INTERN_(kw,child2),this->PonChildren[2]);
  node->field_if_not_unbound(INTERN_(kw,child3),this->PonChildren[3]);
  node->field_if_not_unbound(INTERN_(kw,child4),this->PonChildren[4]);
  node->field_if_not_unbound(INTERN_(kw,child5),this->PonChildren[5]);
  node->field_if_not_unbound(INTERN_(kw,child6),this->PonChildren[6]);
  node->field_if_not_unbound(INTERN_(kw,child7),this->PonChildren[7]);
}

string OctNode_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->_instanceClass()->_classNameAsString();
  ss << " :corner [ " << this->vCorner.getX() << " " << this->vCorner.getY() << " " << this->vCorner.getZ() << "]";
  ss << " :depth " << this->iDepth;
  ss << " :status " << this->iStatus;
  ss << " :nodenum " << this->iNodeNum;
  ss << "#>";
  return ss.str();
}



void Octree_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,type), this->type);
  node->field( INTERN_(kw,maxdepth), this->_iMaxDepth);
  node->field( INTERN_(kw,treepoints), this->iTreePoints);
  node->field( INTERN_(kw,dielectric), this->iDielectric);
  node->field( INTERN_(kw,gridsize), this->dGridSize);
  node->field_if_not_empty( INTERN_(kw,densities), this->PiDensities);
  node->field_if_not_empty( INTERN_(kw,charges), this->_PfCharges);
  node->field_if_not_empty( INTERN_(kw,halfEdges), this->PdHalfEdges);
  node->field_if_not_empty( INTERN_(kw,halfDiagonals), this->PdHalfDiagonals);
  node->field_if_not_empty( INTERN_(kw,atoms), this->vaAtoms);
  node->field( INTERN_(kw,onHead), this->onHead );
}


double Octree_O::dDistanceSq( Vector3 Pv1, Vector3 Pv2 )
{
  double  x, y, z;
  x = Pv1.getX() - Pv2.getX();
  x = x * x;
  
  y = Pv1.getY() - Pv2.getY();
  y = y * y;
  
  z = Pv1.getZ() - Pv2.getZ();
  z = z * z;
  
  return(x + y + z);
}

double Octree_O::dDistance( Vector3 Pv1, Vector3 Pv2 )
{
double  x, y, z;
 
 x = Pv1.getX() - Pv2.getX();
 x = x * x;
 
 y = Pv1.getY() - Pv2.getY();
 y = y * y;
 
 z = Pv1.getZ() - Pv2.getZ();
 z = z * z;
 
 return(sqrt(x + y + z));
}

//static void Octree_O::dumpNode( OctNode_sp PonNode )
void Octree_O::dumpNode( OctNode_sp PonNode)
{
  int	i;
  fprintf(stderr, "@");
  for (i=0; i<PonNode->iDepth; i++)
    fprintf(stderr, " ");
  fprintf(stderr, "%d %d %f %f %f\n", 
          PonNode->iDepth, PonNode->iStatus,
          PonNode->vCorner.getX(), PonNode->vCorner.getY(), PonNode->vCorner.getZ());
  if ( PonNode->iDepth == this->_iMaxDepth  &&  PonNode->iStatus == OCT_PARTIAL) {
    fprintf(stderr, " bad node\n");
#if (!defined WIN32)
		//kill(0,5);
#endif
  }
//	if ( PonNode->iStatus == OCT_PARTIAL  &&  PonNode->PonChildren == NULL ) {
  for (i=0; i<8; i++){
    if ( PonNode->iStatus == OCT_PARTIAL  &&  PonNode->PonChildren[i].nilp() ) {
      fprintf(stderr, "partial without children\n");
#if (!defined WIN32)
		//kill(0,5);
#endif
    }
//	if ( PonNode->PonChildren != NULL ) {
    if ( PonNode->PonChildren[i].notnilp() ) {
      if ( PonNode->iDepth >= this->_iMaxDepth ) {
        fprintf(stderr, "children too deep\n");
#if (!defined WIN32)
			//kill(0,5);
#endif
      }
	//	for (i=0; i<8; i++) 
      dumpNode( PonNode->PonChildren[i]);
    }
  }
}

void Octree_O::dumpOctree()// Octree_sp octTree )
{
  dumpNode( this->onHead);
}

/*************************************************************************
 *************************************************************************/

//static OctNode_sp Octree_O::PonMakeChildren( PonNode_sp ponNode, int iDepth, int iStatus )
void Octree_O::PonMakeChildren(OctNode_sp ponNode, int iDepth, int iStatus, size_t& iNodeCount )
{
//OctNode_sp      PonChildren[8];
  int	       	i;
  double		dHalfEdge = this->PdHalfEdges[iDepth-1];

  for (i=0; i<8; i++) {
    ponNode->PonChildren[i] = OctNode_O::create();
    ponNode->PonChildren[i]->iStatus = iStatus;
    ponNode->PonChildren[i]->iDepth = iDepth;
    ponNode->PonChildren[i]->iAtoms = 0;
    ponNode->PonChildren[i]->vCorner = ponNode->vCorner;
    ponNode->PonChildren[i]->iNodeNum = iNodeCount++;
    
  }
    
	/*
	 *  Make nodes & set simple stuff by initializing 1st node
	 *	and copying it.
	 */
	//GC_ALLOCATE(OctNode_O, PonChildren); //MALLOC( PonChildren, Octnode_sp *, 8 * sizeof(OctNode_sp) );
	//memset(PonChildren, 0, sizeof(OctNode_sp));
	//ponNode->PonChildren->iStatus = iStatus;
	//ponNode->PonChildren->iDepth = iDepth;
	//ponNode->PonChildren->ponNode->PonChildren = NULL;
	//ponNode->PonChildren->PaAtomList = NULL;
	//ponNode->PonChildren->PfCharges = NULL;
	//ponNode->PonChildren[0]->iAtoms = 0;
	//ponNode->PonChildren[0]->vCorner = vCorner;
	//ponNode->PonChildren[0]->iNodeNum = iNodeCount++;
//  ponNode->PonChildren[0]->iStatus = iStatus;
//  ponNode->PonChildren[0]->iDepth = iDepth;
	//ponNode->PonChildren[0]->PonChildren = NULL;
	//ponNode->PonChildren[0]->PaAtomList = NULL;
	//ponNode->PonChildren[0]->PfCharges = NULL;
//  ponNode->PonChildren[0]->iAtoms = 0;
//  ponNode->PonChildren[0]->vCorner = ponNode->vCorner;
//  ponNode->PonChildren[0]->iNodeNum = iNodeCount++;

 // for (i=1; i<8; i++) {

 //   ponNode->PonChildren[i] = ponNode->PonChildren[0];
//		ponNode->PonChildren[i].iNodeNum = iNodeCount++;
//    ponNode->PonChildren[i]->iNodeNum = iNodeCount++;
//  }

	/* 
	 *  Differentiate the non-0-corner children by giving them
	 *	seperate corners of the octant. The order is, do
	 *	Z-stacks along the Y-axis, forming slices on the X-axis.
	 *	(This corresponds to the traversal with inner loop of Z,
	 *	middle loop of Y, and outer loop of X). This order is
	 *	important later for charge redistribution when subdividing
	 *	an existing node to delete a new atom's volume from the tree.
	 */
//	ponNode->PonChildren[1].vCorner.getZ() += dHalfEdge;		/*      0 + Z 	   */
//	ponNode->PonChildren[2].vCorner.getY() += dHalfEdge;		/*      0 + Y 	   */
//
//	ponNode->PonChildren[3].vCorner = ponNode->PonChildren[2].vCorner;
//	ponNode->PonChildren[3].vCorner.getZ() += dHalfEdge;		/*    0 + Y + Z    */
//
//	ponNode->PonChildren[4].vCorner.getX() += dHalfEdge;		/*      0 + X 	   */
//
//	ponNode->PonChildren[5].vCorner = ponNode->PonChildren[4].vCorner;
//	ponNode->PonChildren[5].vCorner.getZ() += dHalfEdge;		/*    0 + X + Z	   */
//	
//	ponNode->PonChildren[6].vCorner = ponNode->PonChildren[4].vCorner;
//	ponNode->PonChildren[6].vCorner.getY() += dHalfEdge;		/*    0 + X + Y	   */
//
//	ponNode->PonChildren[7].vCorner = ponNode->PonChildren[6].vCorner;
//	ponNode->PonChildren[7].vCorner.getZ() += dHalfEdge;		/*  0 + X + Y + Z  */
//  printf("Before addToZ(%lf) PonChildren[0]->vCorner:x  %5.2f  y  %5.2f z %5.2f\n", dHalfEdge, 
//         ponNode->PonChildren[0]->vCorner.getX(), 
//         ponNode->PonChildren[0]->vCorner.getY(), 
//         ponNode->PonChildren[0]->vCorner.getZ());
  ponNode->PonChildren[1]->vCorner.addToZ(dHalfEdge);		/*      0 + Z 	   */
//  printf("After addToZ(%lf) PonChildren[0]->vCorner:x  %5.2f  y  %5.2f z %5.2f\n", dHalfEdge,
//         ponNode->PonChildren[0]->vCorner.getX(), 
//         ponNode->PonChildren[0]->vCorner.getY(), 
//         ponNode->PonChildren[0]->vCorner.getZ());
//  ponNode->PonChildren[2]->vCorner.getY() += dHalfEdge;		/*      0 + Y 	   */
  ponNode->PonChildren[2]->vCorner.addToY(dHalfEdge);		/*      0 + Y 	   */
  
  ponNode->PonChildren[3]->vCorner = ponNode->PonChildren[2]->vCorner;
//  ponNode->PonChildren[3]->vCorner.getZ() += dHalfEdge;		/*    0 + Y + Z    */
  ponNode->PonChildren[3]->vCorner.addToZ(dHalfEdge);		/*    0 + Y + Z    */
  
//  ponNode->PonChildren[4]->vCorner.getX() += dHalfEdge;		/*      0 + X 	   */
  ponNode->PonChildren[4]->vCorner.addToX(dHalfEdge);		/*      0 + X 	   */
  
  ponNode->PonChildren[5]->vCorner = ponNode->PonChildren[4]->vCorner;
//  ponNode->PonChildren[5]->vCorner.getZ() += dHalfEdge;		/*    0 + X + Z	   */
  ponNode->PonChildren[5]->vCorner.addToZ(dHalfEdge);		/*    0 + X + Z	   */
  
  ponNode->PonChildren[6]->vCorner = ponNode->PonChildren[4]->vCorner;
//  ponNode->PonChildren[6]->vCorner.getY() += dHalfEdge;		/*    0 + X + Y	   */
  ponNode->PonChildren[6]->vCorner.addToY(dHalfEdge);		/*    0 + X + Y	   */
  
  ponNode->PonChildren[7]->vCorner = ponNode->PonChildren[6]->vCorner;
//  ponNode->PonChildren[7]->vCorner.getZ() += dHalfEdge;		/*  0 + X + Y + Z  */
  ponNode->PonChildren[7]->vCorner.addToZ(dHalfEdge);		/*  0 + X + Y + Z  */
  
//	return( PonChildren );
}

//static int Octree_O::dFinalCheck( OctNode_sp PonNode, int iAtoms, Atom_sp PaAtomList )
int Octree_O::dFinalCheck( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp> PaAtomList, double dShellRadius )
{
  int	i;
  Atom_sp PaAtom;
#ifdef OCTDEBUG
  depth[PonNode->iDepth]++;
#endif
	/*
	 *  Only need to check corner point. Must be in at least
	 *	one outer radius and in no inner radius.
	 */
  PonNode->iStatus = OCT_EXCLUDED;
  for (i=0; i<iAtoms; i++) {
    PaAtom = PaAtomList[i];
    double	d;
		//d = dDistance( PonNode->vCorner, vAtomPosition( PaAtom ) );
    d = dDistance( PonNode->vCorner, PaAtom->getPosition() );
    if ( d < dShellRadius )
      PonNode->iStatus = OCT_INCLUDED;
		//if ( d < dAtomTemp(PaAtom)) {
    if ( d < PaAtom->dAtomTemp) {
#ifdef OCTDEBUG
      insex[oct->depth]++;
#endif
      PonNode->iStatus = OCT_EXCLUDED;
      break;
    }
  }
#ifdef OCTDEBUG
  if (PonNode->iStatus == OCT_INCLUDED)
    included[oct->depth]++;
  else
    outside[oct->depth]++;
#endif
  return(PonNode->iStatus);
}

/*
 *  DestroyOctant() - destroys any type of octree
 */
//static void Octree_O::DestroyOctant( OctNode_sp PonNode, int iStatus )
void Octree_O::DestroyOctant( OctNode_sp PonNode, int iStatus )
{
     printf("%s:%d:%s vPoint -> %d, %d, %d\n", __FILE__, __LINE__, __FUNCTION__, PonNode->iStatus, PonNode->iNodeNum, PonNode->iDepth);

#ifdef DBG2
  fprintf(stderr, "@@@ destroy node 0x%x level %d type %d -> %d", 
          PonNode, PonNode->iLevel,  PonNode->iStatus, iStatus);
  if (PonNode->PaAtomList.notnilp()) //if ( PonNode->PaAtomList != NULL )
    fprintf(stderr, " atomlist 0x%x", PonNode->PaAtomList);
  else
    fprintf(stderr, " atomlist null");
  if (PonNode->Ponchildren.notnilp()) //if ( PonNode->PonChildren != NULL )
    fprintf(stderr, " children 0x%x\n",  PonNode->PonChildren);
  else
    fprintf(stderr, " children null\n");
#endif
  
  if (PonNode->PaAtomList.size()>0) //if ( PonNode->PaAtomList != NULL )
    PonNode->PaAtomList.clear();//	FREE( PonNode->PaAtomList );

//    if(PonNode->PonChildren[j].notnilp()){//if ( PonNode->PonChildren != NULL ) {
  if(PonNode->PonChildren[0]){//if ( PonNode->PonChildren != NULL ) {
      
      int i;
      
		/*
		 *  Destroy children. Recursive destroys completely
		 *	delete nodes, so mark OCT_UNKNOWN to avoid
		 *	confusion..
		 */
      for (i=0; i<8; i++){
        printf("%s:%d:%s vPoint -> %d, %d, %d %d\n", __FILE__, __LINE__, __FUNCTION__, PonNode->iStatus, PonNode->iNodeNum, PonNode->iDepth, i);
        DestroyOctant( PonNode->PonChildren[i], OCT_UNKNOWN );
      }
      for (i=0; i<8; i++){
        PonNode->PonChildren[i]=_Nil<core::T_O>();//	FREE( PonNode->PonChildren );
      }
    }
  PonNode->iStatus = iStatus;
  
}
//static int Octree_O::iBuildShellOctant( OctNode_sp PonNode, int iAtoms, Atom_sp *PaAtomList )
int Octree_O::iBuildShellOctant( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp> PaAtomList, double dShellRadius )
{
  int	         	i, iPartialIn, iPartialOut, iIncluded, iExcluded, iNewAtoms;
  Atom_sp	        PaAtom;
  gctools::Vec0<Atom_sp>  PaNewAtoms;
  double	         	d, dHalfEdge, dHalfDiagonal;
  Vector3	        	vCenter;
 // OctNode_sp       	PonChildren[8];
  
  
	//PonNode->PaAtomList = NULL;
  PonNode->iAtoms = 0;
  
	/*
	 *  End recursion if the limit of resolution (octree max depth)
	 *	has been reached.
	 */

  if ( PonNode->iDepth == this->_iMaxDepth ) 
    return( dFinalCheck( PonNode, iAtoms, PaAtomList, dShellRadius ) );
  
	/*
	 *  At less than maximum depth. Evaluate whole box with respect 
	 *	to atoms by considering overlaps with single atoms' inner 
	 *	and outer radii. Use center of box +- half diagonal for test.
	 *	Optimally, whole box can be included or excluded. Otherwise,
	 *	it will be divided into 8 and the parts evaluated.
	 */
  
  dHalfEdge = this->PdHalfEdges[PonNode->iDepth];
  vCenter = PonNode->vCorner;
  vCenter.getX() += dHalfEdge;
  vCenter.getY() += dHalfEdge;
  vCenter.getZ() += dHalfEdge;
  dHalfDiagonal = this->PdHalfDiagonals[PonNode->iDepth];
  
	/*
	 *  Get memory for list of atoms overlapping this cube, which 
	 *	is necessarily <= size of parent's.
	 */
  
	//GC_ALLOCATE(Atom_O, PaNewAtoms); // MALLOC( PaNewAtoms, ATOM *, iAtoms * sizeof(ATOM) );
	//memset( PaNewAtoms, 0, iAtoms * sizeof(ATOM) );
  //PaNewAtoms.resize(iAtoms, _Nil<Atom_O>());
  //PonNode->PaAtomList = PaNewAtoms;
  iNewAtoms = 0;
#ifdef DBG2
  fprintf(stderr, "@@@ node 0x%x atoms 0x%x\n", PonNode, PaNewAtoms);
#endif
  
	/*
	 *  Consider all parent's atoms in relation to this box.
	 */
  iIncluded = 0;
  iPartialIn = 0;
  iPartialOut = 0;
  for (i=0; i<iAtoms; i++) {
    PaAtom = PaAtomList[i];
    
		//d = dDistance( &vCenter, &vAtomPosition( PaAtom ) );
    d = dDistance( vCenter, PaAtom->getPosition() );
//    printf("dAtomTem %5.2f\n", PaAtom->dAtomTemp  );
//    printf("d + dHalfDiagonal %5.2f\n", d + dHalfDiagonal );
//    printf("dShellRadius %5.2f\n", dShellRadius );
//    printf("d - dHalfDiagonal %5.2f\n", d - dHalfDiagonal );
    
		/*
		 *  Consider inner radius first, since complete
		 *	inclusion in it means the box is excluded
		 *	and the outer radius doesn't need to be looked at.
		 */
		//if ( d + dHalfDiagonal < dAtomTemp( PaAtom ) ) {
    if ( d + dHalfDiagonal < PaAtom->dAtomTemp ) {
			/*
			 *  Complete inclusion in inner radius of
			 *	this atom: no need to ever look at 
			 *	this box again.
			 */
#ifdef OCTDEBUG
      insex[oct->depth]++;
      depth[oct->depth]++;
#endif
#ifdef DBG2
      fprintf(stderr, "@@@ node 0x%x free atoms 0x%x\n", 
              PonNode, PonNode->PaAtomList);
#endif
			//FREE( PaNewAtoms );
			//PonNode->PaAtomList = NULL;
      PonNode->iStatus = OCT_EXCLUDED;
      return(OCT_EXCLUDED);
    } 
    
		/*
		 *  Consider outer radius
		 */
    if ( d - dHalfDiagonal < dShellRadius ) { 
			/* 
			 *  Box is at least partially inside atom's outer 
			 *	radius: add to new list. Even if box
			 *	is completely included in an outer radius, 
			 *	later discovery of partial overlap with
			 *	some atom's inner radius will force deeper
			 *	evaluation, so the list must be complete.
			 */
      PaNewAtoms.push_back(PaAtom);
      PonNode->PaAtomList.push_back(PaAtom);
      iNewAtoms++;
//      PaNewAtoms[iNewAtoms++] = PaAtom;
      
      if ( d + dHalfDiagonal < dShellRadius ){
        iIncluded++; /* completely inside outer */
      }
      else{
        iPartialOut++; /* partially inside outer */
      }
			//if ( d - dHalfDiagonal < dAtomTemp( PaAtom ) ) {
      if ( d - dHalfDiagonal < PaAtom->dAtomTemp ) {

				/* 
				 *  Partial inclusion in inner radius of
				 *	this atom.
				 */
        iPartialIn++;
      }
    }
  }
 	/*
	 *  Box is not completely inside some inner radius. If it
	 *	is partially inside one or more inner radii, it will
	 *	be evaluated deeper. But first, if box is fully included 
	 *	in target zone or completely outside outer radius, finish 
	 *	it off.
	 */
  if (!iPartialIn) {
    if (iIncluded) {
			/*
		 	 *  Box is completely in 1 or more outer radii and
		 	 *	does not overlap any inner radii.
		 	 */
#ifdef OCTDEBUG
      included[oct->depth]++;
      depth[oct->depth]++;
#endif
      PonNode->iStatus = OCT_INCLUDED;
      PonNode->iAtoms = iNewAtoms;
      return(OCT_INCLUDED);
    } 
    if ( !iPartialOut ) {
			/*
			 *  Box is completely outside outer radii
			 */
#ifdef OCTDEBUG
      outside[oct->depth]++;
      depth[oct->depth]++;
#endif
#ifdef DBG2
      fprintf(stderr, "@@@ node 0x%x free atoms 0x%x\n", 
              PonNode, PonNode->PaAtomList);
#endif
      PonNode->iStatus = OCT_EXCLUDED;
			//FREE( PaNewAtoms );
			//PonNode->PaAtomList = NULL;
      return(OCT_EXCLUDED);
    }
  }
  PonNode->iAtoms = iNewAtoms;

	/*
	 *  Partial inclusion in target zone (part(s) of box overlap
	 *	inner radius or outer).
	 *  Create sub-boxes for recursion.
	 */
//	PonChildren = PonMakeChildren( &PonNode->vCorner, 
//					PonNode->iDepth + 1, OCT_UNKNOWN );
  PonMakeChildren(PonNode,PonNode->iDepth + 1, OCT_UNKNOWN, this->iNodeCount );
  
	/*
	 *  Recurse on sub-boxes, keeping track of their status.
	 */
  iIncluded = 0;
  iExcluded = 0;
  for (i=0; i<8; i++) {
    switch (this->iBuildShellOctant( PonNode->PonChildren[i], 
                                     iNewAtoms, PaNewAtoms, dShellRadius ) ) {
    case OCT_INCLUDED:
        iIncluded++;
        break;
    case OCT_EXCLUDED:
        iExcluded++;
        break;
    case OCT_PARTIAL:
        break;
    default:
        SIMPLE_ERROR(BF("bad type\n"));
//        fprintf(stderr, "bad type\n");
//        exit(1);
    }
  }
//    printf("iIncluded %d\n", iIncluded );
//    printf("iExcluded %d\n", iExcluded );

	/*
	 *  See if whole box is inside or outside zone of interest
	 *	as a sum of subbox-singleatom considerations.
	 */
  
  if ( iExcluded == 8 ) {
    
		/*
		 *  Whole box is within vdws of multiple atoms.
		 */
    
#ifdef OCTDEBUG
    multex[oct->depth]++;
#endif
 //   DestroyOctant( PonNode, OCT_EXCLUDED );
    PonNode->iStatus = OCT_EXCLUDED;
    return( OCT_EXCLUDED );
  } 
  if ( iIncluded == 8 ) {
    
		/*
		 *  Whole box included via outer radii of multiple atoms.
		 */
#ifdef OCTDEBUG
    multwhole[oct->depth]++;
#endif
		//FREE( PonChildren );
		//PonNode->PonChildren = NULL;
    for (i=0; i<8; i++) {
      PonNode->PonChildren[i] = _Nil<core::T_O>();
    }
    PonNode->iStatus = OCT_INCLUDED;
    return(OCT_INCLUDED);
  }
	/*
	 *  Remaining case: some children included, others not.
	 */
  this->fVolume += iIncluded * dHalfEdge * dHalfEdge * dHalfEdge;
//  iTreeGridPoints += iIncluded * PiDensities[PonNode->iDepth+1]; 
  this->iTreePoints += iIncluded * this->PiDensities[PonNode->iDepth+1]; 
  PonNode->iStatus = OCT_PARTIAL;
  return(OCT_PARTIAL);
}

/*
 *  For an interior octant, the point is to keep atom lists at
 *	the bottom nodes of the tree.
 */
//static int Octree_O::iBuildInteriorOctant( OctNode_sp PonNode, int iAtoms, Atom_sp *PaAtomList )
int Octree_O::iBuildInteriorOctant( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp> PaAtomList)
{
  int		         i, iIncluded, iExcluded; 
  Atom_sp		 PaAtom;
  gctools::Vec0<Atom_sp> PaNewAtoms;
  int		         iNewAtoms;
  double		 d, dCenterRadiusSq, dHalfEdge;
  Vector3                vCenter;
  OctNode_sp	         PonChildren[8];
    
#ifdef OCTDEBUG
  depth[PonNode->iDepth]++;
#endif
  
	/*
	 *  Evaluate whole box with respect to atoms by considering
	 *	overlaps with single atoms' inner and outer radii.
	 *	Use center of box + half length of diagonal for test.
	 */
  
  dHalfEdge = this->PdHalfEdges[PonNode->iDepth];
  vCenter = PonNode->vCorner;
  vCenter.getX() += dHalfEdge;
  vCenter.getY() += dHalfEdge;
  vCenter.getZ() += dHalfEdge;
  dCenterRadiusSq = this->PdHalfDiagonals[PonNode->iDepth];
  dCenterRadiusSq += 3; /* TODO - arbitrary.. */
	/* dCenterRadiusSq *= dCenterRadiusSq; */
  
	/*
	 *  Get memory for list of atoms inside sphere enclosing this cube, which 
	 *	must be <= size of parent's.
	 */
	//GC_ALLOCATE(Atom_O, PaNewAtoms); // MALLOC( PaNewAtoms, ATOM *, iAtoms * sizeof(ATOM) );
	//memset( PaNewAtoms, 0, iAtoms * sizeof(ATOM) );
  //PonNode->PaAtomList = PaNewAtoms;
#ifdef DBG2
  fprintf(stderr, "@@@ inode 0x%x atoms 0x%x\n", PonNode, PonNode->PaAtomList);
#endif
  iNewAtoms = 0;
  
	/*
	 *  Consider all parent's atoms in relation to this box.
	 */
  iIncluded = 0;
  for (i=0; i<iAtoms; i++) {
    PaAtom = PaAtomList[i];
		//d = dDistance( &vCenter, &vAtomPosition( *PaAtom ) );
                //d = dDistance( vCenter, PaAtom->getPosition() );
    if ( d < dCenterRadiusSq ) {
			/* 
			 *  It's w/in larger sphere
			 */
      PaNewAtoms.push_back(PaAtom);
      PonNode->PaAtomList.push_back(PaAtom);
      iNewAtoms++;

//      PaNewAtoms[iNewAtoms++] = PaAtom;
    }
  }
  if (!iNewAtoms) {
		/*
		 *  No atoms w/in box-enclosing sphere
		 */
#ifdef OCTDEBUG
    outside[oct->depth]++;
    depth[oct->depth]++;
#endif
		//FREE( PonNode->PaAtomList );
		// PonNode->PaAtomList = NULL;
    PonNode->iStatus = OCT_EXCLUDED;
    return(OCT_EXCLUDED);
  }
  
  PonNode->iAtoms = iNewAtoms;
  
	/*
	 *  End recursion if the limit of resolution (octree max depth)
	 *	has been reached or if the density is adequate for 
	 *	linear search.
	 */
  if ( PonNode->iDepth == this->_iMaxDepth || ( iNewAtoms && iNewAtoms < 10 )) {
//    iTreeGridPoints += PiDensities[PonNode->iDepth]; 
    this->iTreePoints += this->PiDensities[PonNode->iDepth]; 
    PonNode->iStatus = OCT_INCLUDED;
    return(OCT_INCLUDED);
  }
  
	/*
	 *  Create sub-boxes for recursion.
	 */
  PonMakeChildren( PonNode, PonNode->iDepth + 1, OCT_UNKNOWN, this->iNodeCount );
  //PonNode->PonChildren = PonChildren;
  
	/*
	 *  recurse on sub-boxes
	 */
  iIncluded = 0;
  iExcluded = 0;
  for (i=0; i<8; i++) {
    switch (iBuildInteriorOctant( PonChildren[i], 
                                  iNewAtoms, PaNewAtoms ) ) {
      
    case OCT_INCLUDED:
        iIncluded++;	/* for debug tracking */
        break;
    case OCT_EXCLUDED:
        iExcluded++;
        break;
    default:
        SIMPLE_ERROR(BF("bad type\n"));
//        fprintf(stderr, "bad type\n");
//        exit(1);
    }
  }
  
	/*
	 *  
	 */
  if ( iExcluded == 8 ) {
		//FREE( PonNode->PonChildren );
    for (i=0; i<8; i++) {
 		//PonNode->PonChildren = NULL;
      PonNode->PonChildren[i] = _Nil<core::T_O>();
    }
  } else {
		//FREE( PonNode->PaAtomList );
#ifdef DBG2
    fprintf(stderr, "@@@ node 0x%x free atoms 0x%x\n", PonNode, PonNode->PaAtomList);
#endif
		//PonNode->PaAtomList = NULL;
   // PonNode->PaAtomList = _Nil<core::T_O>();
  }
  
  PonNode->iStatus = OCT_INCLUDED;
  return(OCT_INCLUDED);
}

SYMBOL_EXPORT_SC_(KeywordPkg,solvent);
SYMBOL_EXPORT_SC_(ChemPkg,lookup_atom_properties_radius);
//CL_DEFMETHOD Octree_sp Octree_O::octOctTreeCreate(Aggregate_sp uUnit, int iType, double dGridSpace, double dAddExtent, double dShellExtent,
CL_DEFMETHOD void Octree_O::octOctTreeCreate(Aggregate_sp uUnit, OctreeType type, double dGridSpace, double dAddExtent, double dShellExtent, FFNonbondDb_sp nonbondDb,  int iIncludeSolvent, bool bVerbose)
{
 // Octree_sp		  octTree;
  Vector3                 vMinCorner, vMaxCorner, vAtom;
  gctools::Vec0<Atom_sp>  vaAtoms;
  Atom_sp		  aAtom;
  //gctools::Vec0<Atom_sp>  PaAtoms;
  Loop		          lAtoms, lRes;
  Residue_sp		  rRes;
  size_t                   imd;
  int	          	i, j, iAtoms, iDefaultedRadius, iBuild;
  double		dMaxRadius, dCharge, dShellRadius;
  double		dTx, dTy, dTz, dTmax, dTmp, volumePercentage;

  ASSERT(uUnit && uUnit.notnilp());

	/*
	 *  Set globals.
	 */
//  dGridSize = dGridSpace;
  dShellRadius = dShellExtent;	/* shell: clearance beyond inner shell*/

	/*
	 *  Create the octree "object" and initialize
	 */
  //GC_ALLOCATE(Octree_O, octTree); // MALLOC( octTree, OCTREE, sizeof(OCTREEt) );
 // octTree->iType = iType;
  //octTree->dGridSize = dGridSpace;
  this->type = type;
  this->dGridSize = dGridSpace;
 // octTree->PfCharges = NULL;
  // octTree->PdHalfEdges = NULL;
  // octTree->PdHalfDiagonals = NULL;

	/*
	 *  Make array for atom pointers
	 */
	// vaAtoms = vaVarArrayCreate( sizeof(ATOM) );

	/*
	 *  Fill array of atom pointers according to octree type, initializing
	 *	the atoms' temporary floating-point values w/ their sizes.
	 */

  iDefaultedRadius = 0;
  lRes.loopTopAggregateGoal(uUnit, RESIDUES); // lRes = lLoop( uUnit, RESIDUES );
  switch ( type ) {
  case Shell:
  case InteriorSolute:
      // while ((rRes = (RESIDUE) oNext(&lRes))) {
      while (lRes.advance()) {
        rRes = lRes.getResidue();
        if ( iIncludeSolvent  ||  
             rRes->_Type != kw::_sym_solvent){// cResidueType( rRes ) != RESTYPESOLVENT ) {
          lAtoms.loopTopGoal(rRes,ATOMS); // lAtoms = lLoop( (OBJEKT)rRes, ATOMS );
          while (lAtoms.advance()) { // while ((aAtom = (ATOM) oNext(&lAtoms))) {
            aAtom = lAtoms.getAtom();
            vaAtoms.push_back(aAtom); // VarArrayAdd( vaAtoms, (GENP)&aAtom );
            // iDefaultedRadius += iAtomSetTmpRadius( aAtom );
            core::T_mv result = core::eval::funcall(chem::_sym_lookup_atom_properties_radius,aAtom, nonbondDb);
            aAtom->dAtomTemp = core::clasp_to_double(result);
            core::T_sp defaulted = result.second();
            iDefaultedRadius += defaulted.isTrue() ? 1 : 0;
          }
        }
      }
      break;

  case InteriorSolvent: {
    int nowarning = 1;
    //while ((rRes = (RESIDUE) oNext(&lRes))) {
    while (lRes.advance()){
      rRes = lRes.getResidue();
      if ( rRes->_Type == kw::_sym_solvent ) {
        i = 0;
        lAtoms.loopTopGoal(rRes,ATOMS); //lAtoms = lLoop( (OBJEKT)rRes, ATOMS );
        while (lAtoms.advance()) { //while ((aAtom = (ATOM) oNext(&lAtoms))) {
          if (!i)
            vaAtoms.push_back(aAtom); // VarArrayAdd( vaAtoms, (GENP)&aAtom );
          //iDefaultedRadius += iAtomSetTmpRadius( aAtom );
          core::T_mv result = core::eval::funcall(chem::_sym_lookup_atom_properties_radius,aAtom);
          aAtom->dAtomTemp = core::clasp_to_double(result);
          core::T_sp defaulted = result.second();
          iDefaultedRadius += defaulted.isTrue() ? 1 : 0;
          i++;
        }
        if ( i > 3  &&  nowarning ) {
          printf( "Warning: non-water solvent may lead to steric problems");
          nowarning = 0;
        }
      }
    }
  }
      break;
  default:
      SIMPLE_ERROR(BF("Octree type not implemented\n"));
//      printf("Octree type not implemented\n");
//      exit(1);
  }

	/*
	 *  If no atoms, nothing to build.
	 */
  iAtoms = vaAtoms.size(); // iVarArrayElementCount( vaAtoms );
/*
VP0(("atoms: %d\n", iAtoms));
dTmax = 20.0;
for (dTmp=dTmax, iMaxDepth=0; dTmp>dGridSize; iMaxDepth++, dTmp/=2.0);
octTree->iMaxDepth = iMaxDepth;
for (dTmax=dGridSize, i=0; i<iMaxDepth; i++, dTmax*=2.0);
MALLOC( PdHalfEdges, double *, (iMaxDepth+1) * sizeof(double) );
octTree->PdHalfEdges = PdHalfEdges;
MALLOC( PdHalfDiagonals, double *, (iMaxDepth+1) * sizeof(double) );
octTree->PdHalfDiagonals = PdHalfDiagonals;
return(octTree);
*/

  if ( !iAtoms ) {
    //VarArrayDestroy( &vaAtoms );
    //FREE( octTree );
    //return(NULL);
//    return(_Nil<core::T_O>());
    SIMPLE_ERROR(BF("There are not atoms, iAtoms==0"));
  }

//  octTree->vaAtoms = vaAtoms;
  this->vaAtoms = vaAtoms;

  if ( iDefaultedRadius )
    printf( "Used default radius \n");
    //      ATOM_DEFAULT_RADIUS, iDefaultedRadius );

	/*
	 *  Find bounding box of atom centers and max radius, while
	 *	adding the additional extents to the temporary atom radii.
	 *	For purposes of finding max radius, initialize 1st atom
	 *	then loop over others, comparing.
	 */
  // PaAtoms = PVAI( vaAtoms, ATOM, 0 );
  dCharge = vaAtoms[0]->charge; // dAtomCharge(  *PaAtoms );
  dMaxRadius = vaAtoms[0]->dAtomTemp; // dAtomTemp( *PaAtoms );
  //vMinCorner = vMaxCorner = vAtomPosition( vaAtoms[0] ); // *PaAtoms );
  vMinCorner = vMaxCorner = vaAtoms[0]->getPosition(); // *PaAtoms );
  vaAtoms[0]->dAtomTemp += dAddExtent; // AtomTempDoubleIncrement( *PaAtoms, dAddExtent );

  for (i=1; i<iAtoms; i++ ) { // ++PaAtoms ) {
    dCharge += vaAtoms[i]->charge; // dAtomCharge( *PaAtoms );
    dMaxRadius = MAX( dMaxRadius, vaAtoms[i]->dAtomTemp); // MAX( dMaxRadius, dAtomTemp( *PaAtoms ) );
    vaAtoms[i]->dAtomTemp += dAddExtent; // AtomTempDoubleIncrement( *PaAtoms, dAddExtent );
    vAtom = vaAtoms[i]->getPosition(); // vAtom = vAtomPosition( *PaAtoms );
    if (vAtom.getX() < vMinCorner.getX())
      vMinCorner.getX() = vAtom.getX();
    else if (vAtom.getX() > vMaxCorner.getX())
      vMaxCorner.getX() = vAtom.getX();
    if (vAtom.getY() < vMinCorner.getY())
      vMinCorner.getY() = vAtom.getY();
    else if (vAtom.getY() > vMaxCorner.getY())
      vMaxCorner.getY() = vAtom.getY();
    if (vAtom.getZ() < vMinCorner.getZ())
      vMinCorner.getZ() = vAtom.getZ();
    else if (vAtom.getZ() > vMaxCorner.getZ())
      vMaxCorner.getZ() = vAtom.getZ();
  }

	/*
	 *  Expand bounding box if neccessary.
	 */
  switch ( type ) {
  case Shell:
		/*
		 *  Offset shell by max solute atom radius
		 */
      dShellRadius += dMaxRadius + dAddExtent;
      
		/*
		 *  Enclose atom centers with shell
		 */
      vMinCorner.getX() -= dShellRadius;
      vMinCorner.getY() -= dShellRadius;
      vMinCorner.getZ() -= dShellRadius;
      vMaxCorner.getX() += dShellRadius;
      vMaxCorner.getY() += dShellRadius;
      vMaxCorner.getZ() += dShellRadius;
      break;
      
  case InteriorSolute:
  case InteriorSolvent:
		/*
		 *  Enclose atom centers with radius of largest atom
		 */
      vMinCorner.getX() -= dMaxRadius;
      vMinCorner.getY() -= dMaxRadius;
      vMinCorner.getZ() -= dMaxRadius;
      vMaxCorner.getX() += dMaxRadius;
      vMaxCorner.getY() += dMaxRadius;
      vMaxCorner.getZ() += dMaxRadius;
      break;
  default:
      SIMPLE_ERROR(BF("Octtree type %d is not implemented") % type);
  }

  if (bVerbose) {
    core::write_bf_stream(BF("Total solute charge:  %5.2f  Max atom radius:  %5.2f\n") % dCharge % dMaxRadius );
    if ( type == Shell )
      core::write_bf_stream(BF("Grid extends from solute vdw + %.2f  to  %.2f\n") % dAddExtent % dShellRadius );
    core::write_bf_stream(BF("Box:\n" ));
    core::write_bf_stream(BF("   enclosing:  %5.2f %5.2f %5.2f   %5.2f %5.2f %5.2f\n") %  vMinCorner.getX() % vMinCorner.getY() % vMinCorner.getZ()
              % vMaxCorner.getX() % vMaxCorner.getY() % vMaxCorner.getZ());
  }

	/*  
	 *  Find longest edge of box
	 */
  dTx = vMaxCorner.getX() - vMinCorner.getX();
  dTy = vMaxCorner.getY() - vMinCorner.getY();
  dTz = vMaxCorner.getZ() - vMinCorner.getZ();

  dTmax = MAX( dTx, dTy );
  dTmax = MAX( dTmax, dTz );

	/*  
	 *  Determine octree descent depth required to obtain 
	 *	requested grid resolution by dividing longest edge in
	 *	half until it is less than dGridSize.
	 */
  for (dTmp=dTmax, imd=0; dTmp>this->dGridSize; imd++, dTmp/=2.0);
 // octTree->iMaxDepth = iMaxDepth;
  this->_iMaxDepth = imd;

	/*  
	 *  Scale box up to modulo grid size 
	 */
  for (dTmax=this->dGridSize, i=0; i<imd; i++, dTmax*=2.0);

	/*
	 *  Just for fun, calc the new max corner
	 */
  vMaxCorner.getX() += dTmax - dTx;
  vMaxCorner.getY() += dTmax - dTy;
  vMaxCorner.getZ() += dTmax - dTz;

	/*  
	 *  Set up arrays (indexed by depth) of edge size & half-diagonals 
	 *	to save recalculation
	 */
//  PdHalfEdges.resize(imd+1,0.0); // MALLOC( PdHalfEdges, double *, (iMaxDepth+1) * sizeof(double) );
//  octTree->PdHalfEdges = PdHalfEdges;
  this->PdHalfEdges.resize(imd+1,0.0);
//  PdHalfDiagonals.resize(imd+1,0.0); // MALLOC( PdHalfDiagonals, double *, (iMaxDepth+1) * sizeof(double) );
//  octTree->PdHalfDiagonals = PdHalfDiagonals;
  this->PdHalfDiagonals.resize(imd+1,0.0);

  for (dTmp=dTmax/2.0, i=0; i<=imd; dTmp/=2.0, i++) {
    this->PdHalfEdges[i] = dTmp;
		/* diagonal = sqrt( 3 * side^2 ) */
    this->PdHalfDiagonals[i] = sqrt( 3.0 * (dTmp * dTmp) );
    if(bVerbose){
      core::write_bf_stream(BF("PdHalfEdges:        %5.2f\n") % this->PdHalfEdges[i] );
      core::write_bf_stream(BF("PdHalfDiagonals:     %5.2f Angstrom.\n") % this->PdHalfDiagonals[i] );
    }

  }

	/*  
	 *  Set up array (indexed by depth) of points per box, used
	 *	for rejuggling charge array when splitting boxes.
	 */
 // PiDensities.resize(imd+1,0.0); // MALLOC( PiDensities, int *, (iMaxDepth+1) * sizeof(int) );
//  octTree->PiDensities = PiDensities;
  this->PiDensities.resize(imd+1,0.0);
  for (i=imd, j=1; i>-1; i--, j*=8)
     this->PiDensities[i] = j;
  
  if (bVerbose) {
    core::write_bf_stream(BF("   sized:\t\t\t      %5.2f %5.2f %5.2f\n") % vMaxCorner.getX() % vMaxCorner.getY() % vMaxCorner.getZ());
    core::write_bf_stream(BF("   edge:        %5.2f\n") % dTmax );
    core::write_bf_stream(BF("Resolution:     %5.2f Angstrom.\n") % this->dGridSize );
    core::write_bf_stream(BF("Tree depth: %d\n") % imd);
  }
	/*
	 *  Build head node w/ all atoms in list.
	 */
//  octTree->onHead->iStatus = OCT_UNKNOWN;
//  octTree->onHead->iDepth = 0;
//  octTree->onHead->vCorner = vMinCorner;
  this->onHead = OctNode_O::create();
  this->onHead->iStatus = OCT_UNKNOWN;
  this->onHead->iDepth = 0;
  this->onHead->vCorner = vMinCorner;

	/* 
	 *  The master atom list is part of a vararray under the OCTREE
	 *	so, unlike the seperately allocated atom lists in the
	 *	lower nodes on the tree, the master list is not freed:
	 *	rather, its vararray is destroyed. Reset the ATOM* ptr.
	 */
  //octTree->onHead.PaAtomList = NULL;	
  //PaAtoms = PVAI( vaAtoms, ATOM, 0 );

	/*
	 *  Build octree recursively.
	 */
//  iTreeGridPoints = 0;
  this->fVolume = 0.0;
  this->iNodeCount = 0;
//  time_start = time((time_t *) 0);

  switch ( type ) {
  case Shell:
      iBuild = iBuildShellOctant( this->onHead, iAtoms, vaAtoms, dShellRadius ); //PaAtoms );
      break;
  case InteriorSolute:
  case InteriorSolvent:
      iBuild = iBuildInteriorOctant( this->onHead, iAtoms, vaAtoms ); //PaAtoms );
      break;
  default:
      SIMPLE_ERROR(BF("bad switch\n"));
//      printf(("bad switch\n"));
//      exit(1);
  }
//  octTree->iTreePoints = iTreeGridPoints;	/* global */
//  this->iTreePoints = iTreeGridPoints;	/* global */

//  printf( "grid build: %ld sec\n", 
//        time((time_t *)0) - time_start );
  if (bVerbose) {
   volumePercentage = 100 * this->fVolume / ( dTmax * dTmax * dTmax );
   core::write_bf_stream(BF( "Volume = %5.2f of box, grid points %d\n") % volumePercentage % this->iTreePoints ); //iTreeGridPoints );
  }
#ifdef OCTDEBUG
  printf(stderr, "depth  r_inc   r_inex  r_out  multin  multout \n");
  for (i=0;i<=imd; i++) {
    printf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\n", 
            depth[i], 
            included[i],
            insex[i], 
            outside[i], 
            multwhole[i],
            multex[i]
            );
  }

#endif
	/*
	 *  Reset atom radii.
	 */
  //PaAtoms = PVAI( vaAtoms, ATOM, 0 );
  for (i=0; i<iAtoms; i++)
    vaAtoms[i]->dAtomTemp -= dAddExtent; //AtomTempDoubleIncrement( *PaAtoms, -dAddExtent );

  //return(this->asSmartPtr());
  return;
}

/*************************************************************************

BuildShellOctant: Build octree for shell around a list of atoms. 

	"Included" means in the shell.
	"Excluded" means either inside an atom or beyond the shell.

	Lists of atoms enclosing or intersecting a box with outer radius
	are built and passed to the children to reduce their evaluation 
	time. For the shell octree, these lists are not saved.

 *************************************************************************/





/*************************************************************************
 *************************************************************************/

void Octree_O::OctTreeDestroy()// Octree_sp PoctTree )
{
	/*
	 *  Free top-level things.
	 */
  //if ( PoctTree->vaAtoms )
    //VarArrayDestroy( PoctTree->vaAtoms );
  
//	if ( (*PoctTree)->PfCharges != NULL )
	//if ( (*PoctTree)->PfCharges.notnilp())
		//FREE( (*PoctTree)->PfCharges );
  
	//FREE( (*PoctTree)->PdHalfEdges );
	//FREE( (*PoctTree)->PdHalfDiagonals );
	//FREE( (*PoctTree)->PiDensities );

	/*
	 *  Set global and recursively delete the tree.
	 */
  //size_t imd = PoctTree->_iMaxDepth;
  //DestroyOctant( PoctTree->onHead, OCT_UNKNOWN );
	//FREE( *PoctTree );
	//*PoctTree = NULL;
  this->asSmartPtr()  = _Nil<core::T_O>();
//  PoctTree = _Nil<core::T_O>();
  return;
}


/* @@@ */
/*************************************************************************
 *************************************************************************/

//static void Octree_O::OctNodeInitCharges( OctNode_sp PonNode )
void Octree_O::OctNodeInitCharges( OctNode_sp PonNode, int iDistanceCharge)
{
  int	i, j, k, l, pfccount;
  int 	ct = this->_iMaxDepth - PonNode->iDepth + 1;
  int	iCompCharge;
  Vector3	vPoint;
  Atom_sp	PaAtom;
  double	d;
  
  if ( PonNode->iStatus == OCT_PARTIAL ) {
    for (i=0; i<8; i++)
      this->OctNodeInitCharges( PonNode->PonChildren[i], iDistanceCharge);
    return;
  }

  if ( PonNode->iStatus == OCT_EXCLUDED )
    return;
  
	/*
	 *  Node is included: calculate charges.
	 * 	Set node pointer into master charge array 
	 */
//  PonNode->_PfCharges = PfCharges;
	/* 
	 *  For each point x,y,z in this quadrant, 
	 *	loop over atoms, accumulating charge 
	 */
  vPoint.getX() = PonNode->vCorner.getX();

  pfccount = 0;
  for (i=0; i<ct; i++, vPoint.getX()+=this->dGridSize) {
    vPoint.getY() = PonNode->vCorner.getY();
    for (j=0; j<ct; j++, vPoint.getY()+=this->dGridSize) {
      vPoint.getZ() = PonNode->vCorner.getZ();
      for (k=0; k<ct; k++, vPoint.getZ()+=this->dGridSize) {
			/*
			 *  Got point: loop over atoms, accumulating charge.
			 */
        PonNode->_PfCharges.push_back(0.0);
        iCompCharge = 1;
        for (l=0; l< this->vaAtoms.size(); l++){ //iChargeAtoms; l++) {
          PaAtom = this->vaAtoms[l] ;//PaChargeAtoms[l];
          double 	dX, dY, dZ;
          
				//dX = vPoint.getX() - vAtomPosition(*PaAtom).getX();
          dX = vPoint.getX() - PaAtom->getPosition().getX();
          dX = dX * dX;
				//dY = vPoint.getY() - vAtomPosition(*PaAtom).getY();
          dY = vPoint.getY() - PaAtom->getPosition().getY();
          dY = dY * dY;
				//dZ = vPoint.getZ() - vAtomPosition(*PaAtom).getZ();
          dZ = vPoint.getZ() - PaAtom->getPosition().getZ();
          dZ = dZ * dZ;
          d = dX + dY + dZ;
          if ( iDistanceCharge )
            d = sqrt(d);
          PonNode->_PfCharges[pfccount] += PaAtom->charge / d;
				//if ( d < dAtomTemp(PaAtom) )
          if ( d < PaAtom->dAtomTemp )
            iCompCharge = 0;
        }
//	printf("%s:%d:%s vPoint -> %d, %d, %d, %f\n", __FILE__, __LINE__, __FUNCTION__, PonNode->iStatus, PonNode->iNodeNum, PonNode->iDepth, PonNode->_PfCharges[pfccount]);
        if ( iCompCharge ) {
				/*
				 *  Keep track of max, min charges and 
				 *	their locations.
				 */
          if ( PonNode->_PfCharges[pfccount] > this->fMaxCharge ) {
            this->fMaxCharge =  PonNode->_PfCharges[pfccount];
            this->vMaxCharge = vPoint;
          } else if (  PonNode->_PfCharges[pfccount] < this->fMinCharge ) {
            this->fMinCharge =  PonNode->_PfCharges[pfccount];
            this->vMinCharge = vPoint;
          }
        }
        pfccount++;
      }
    }
  }

  for (l=0; l< PonNode->_PfCharges.size(); l++) {
    this->_PfCharges.push_back(PonNode->_PfCharges[l]);
    this->iNodeNumCharges.push_back(PonNode->iNodeNum);
  }

  
  return;
}
//CL_DEFMETHOD void Octree_O::OctTreeInitCharges( Octree_sp octTree, int iAtomOption, int iDielectric, double dCutDist,
//                         Vector3& vMin, Vector3& vMax )
CL_DEFMETHOD core::T_mv Octree_O::OctTreeInitCharges( /*Octree_sp octTree,*/ int iAtomOption, int iDielectric, double dCutDist)
{
  int	i, iDistanceCharge, iChargeAtoms;
  Atom_sp	PaAtom;
  
  if ( this->type != Shell ) {
    SIMPLE_ERROR(BF("InitCharges: wrong tree type\n"));
//    printf(( "InitCharges: wrong tree type\n" ));
//    exit(1);
  }

  if ( !this->iTreePoints ) {
    SIMPLE_ERROR(BF("InitCharges: no grid (?)\n"));
//    printf(( "InitCharges: no grid (?)\n" ));
//    exit(1);
  }

  printf( "Calculating grid charges\n" );
  time_start = time((time_t *) 0);
  
	/*
	 *  Get the memory for the charges and note the dielectric.
	 */
  
  //this->_PfCharges.resize(this->iTreePoints,0.0); // MALLOC( PfCharges, float *, octTree->iTreePoints * sizeof(float) );
  this->iDielectric = iDielectric;

	/*
	 *  Set up globals for node descent.
	 */
  if ( iAtomOption == AT_OCTREE ) {
    iChargeAtoms = this->vaAtoms.size(); //iVarArrayElementCount( octTree->vaAtoms );
    //PaChargeAtoms = this->vaAtoms; //PVAI( octTree->vaAtoms, ATOM, 0 );
  } else
    SIMPLE_ERROR(BF("iAtomOption != AT_OCTREE"));
    //exit(1);
  this->fMaxCharge = -FLT_MAX;
  this->fMinCharge = FLT_MAX;
  //dGridSize = octTree->dGridSize;

  
  if ( iDielectric == DIEL_R2 ) {
    iDistanceCharge = 0;
  } else
    iDistanceCharge = 1;
  
	/*
	 *  Add new radius to atoms.
	 */
  if ( iDielectric == DIEL_R2 ) { 
		//PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
    for (i=0; i<iChargeAtoms; i++ ){ //, PaAtom++) {
/* TODO - chargeatoms not necc whole set */
      this->vaAtoms[i]->dAtomTemp += dCutDist; //AtomTempDoubleIncrement( *PaAtom, dCutDist );
      this->vaAtoms[i]->dAtomTemp *= this->vaAtoms[i]->dAtomTemp; //	AtomTempDoubleSquare( *PaAtom );
    }
  } else {
		//PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
    for (i=0; i<iChargeAtoms; i++ ){ //, PaAtom++)
      this->vaAtoms[i]->dAtomTemp += dCutDist; // AtomTempDoubleIncrement( *PaAtom, dCutDist );
    }
  }


	/*
	 *  Descend the octree.
	 */
  this->OctNodeInitCharges( this->onHead, iDistanceCharge);
// vMin = vMinCharge;
// vMax = vMaxCharge;

	/*
	 *  Restore atoms' radii.
	 */
  if ( iDielectric == DIEL_R2 ) { 
		//PaAtom = PVaAI( octTree->vaAtoms, ATOM, 0 );
    for (i=0; i<iChargeAtoms; i++ ){ //, PaAtom++) {
      this->vaAtoms[i]->dAtomTemp = sqrt(this->vaAtoms[i]->dAtomTemp); //AtomTempDoubleSquareRoot( *PaAtom );
      this->vaAtoms[i]->dAtomTemp += -dCutDist; //AtomTempDoubleIncrement( *PaAtom, -dCutDist );
    }
  } else {
		//PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
    for (i=0; i<iChargeAtoms; i++) //, PaAtom++)
      this->vaAtoms[i]->dAtomTemp += -dCutDist; //AtomTempDoubleIncrement( *PaAtom, -dCutDist );
  }
  printf( "charges: %ld sec\n", time((time_t *) 0) - time_start);
  return Values(geom::OVector3_O::create(this->vMinCharge),geom::OVector3_O::create(this->vMaxCharge));
}


/*************************************************************************
 *************************************************************************/

//static void Octree_O::OctNodePrintGrid( OctNode_sp *PonNode )
CL_DEFMETHOD void Octree_O::OctNodePrintGrid( OctNode_sp PonNode, int iColor)
{
  Vector3	vPoint;
  int	i, j, k, ct, ccharge;
  float	PfCharge;
  if ( PonNode->iStatus == OCT_PARTIAL ) {
    for (i=0; i<8; i++) 
      this->OctNodePrintGrid( PonNode->PonChildren[i], iColor );
    return;
  }
  if ( PonNode->iStatus == OCT_EXCLUDED )
    return;
  
	/*
	 *  Included - print the box.
	 */

  ccharge = 0;
  ct = this->_iMaxDepth - PonNode->iDepth + 1;
  vPoint.getX() = PonNode->vCorner.getX();
  for (i=0; i<ct; i++, vPoint.getX()+=this->dGridSize) {
    vPoint.getY() = PonNode->vCorner.getY();
    for (j=0; j<ct; j++, vPoint.getY()+=this->dGridSize) {
      vPoint.getZ() = PonNode->vCorner.getZ();
      for (k=0; k<ct; k++, vPoint.getZ()+=this->dGridSize) {
        PfCharge = PonNode->_PfCharges[ccharge];
//        switch ( iColorMethod ) {
//        case COLOR_RANGE:
//            fprintf(fChargeFile, ".color %d\n",
//                   (int)floor( 5 + 60 *
//                               (PfCharge-this->fMinCharge) / 
//                               (this->fMaxCharge-this->fMinCharge) ) );
//            fprintf(fChargeFile, ".dot %f %f %f\n",
//                   vPoint.getX(), vPoint.getY(), vPoint.getZ());
//            break;
//        case COLOR_CUT:
//            if ( PfCharge < -0.1 )
//              fprintf(fChargeFile, ".color yellow\n");
//            else if ( PfCharge > 0.1 )
//              fprintf(fChargeFile, ".color cyan\n");
//            else
//              fprintf(fChargeFile, ".color black\n");
//            fprintf(fChargeFile, ".dot %f %f %f\n",
//                    vPoint.getX(), vPoint.getY(), vPoint.getZ());
//            break;
//        case COLOR_DEPTH:
//            if ( PonNode->iDepth == iMaxDepth )
//              fprintf(fChargeFile, ".color white\n");
//            else if ( PonNode->iDepth % 2 )
//              fprintf(fChargeFile, ".color red\n");
//            else
//              fprintf(fChargeFile, ".color cyan\n");
//            fprintf(fChargeFile, ".dot %f %f %f\n",
//                   vPoint.getX(), vPoint.getY(), vPoint.getZ());
//            break;
//        case COLOR_NONE:
//            fprintf(fChargeFile, "%f   %f %f %f\n",
//                   PfCharge, 
//                   vPoint.getX(), vPoint.getY(), vPoint.getZ());
//            break;
//        default:
//            fprintf(fChargeFile, ".color white\n");
//        }
        switch ( iColor ) {//switch ( iColorMethod ) {
        case COLOR_RANGE:
            core::write_bf_stream(BF(".color %d\n") % 
                      (int)floor( 5 + 60 *
                                  (PfCharge-this->fMinCharge) / 
                                  (this->fMaxCharge-this->fMinCharge) ) );
            core::write_bf_stream(BF(".dot %f %f %f\n") % 
                      vPoint.getX() % vPoint.getY() % vPoint.getZ());
            break;
        case COLOR_CUT:
            if ( PfCharge < -0.1 ){
              core::write_bf_stream(BF(".color yellow\n"));
            } else if ( PfCharge > 0.1 ){
              core::write_bf_stream(BF(".color cyan\n"));
            } else
              core::write_bf_stream(BF(".color black\n"));
            core::write_bf_stream(BF(".dot %f %f %f\n") % 
                      vPoint.getX() % vPoint.getY() % vPoint.getZ());
            break;
        case COLOR_DEPTH:
            if ( PonNode->iDepth == this->_iMaxDepth ){
              core::write_bf_stream(BF(".color white\n"));
            }else if ( PonNode->iDepth % 2 ){
              core::write_bf_stream(BF(".color red\n"));
            } else
              core::write_bf_stream(BF(".color cyan\n"));
            core::write_bf_stream(BF(".dot %f %f %f\n") % 
                      vPoint.getX() % vPoint.getY() % vPoint.getZ());
            break;
        case COLOR_NONE:
            core::write_bf_stream(BF("%f   %f %f %f\n") % 
                      PfCharge %  
                      vPoint.getX() %  vPoint.getY() %  vPoint.getZ());
            break;
        default:
            core::write_bf_stream(BF(".color white\n"));
        }
        ccharge++ ; //PfCharge++;
      }
    }
  }
}

//void Octree_O::OctTreePrintGrid( Octree_sp octTree, core::T_sp stream, int iColor ) //char *sFileName, int iColor )
//{
//  BFORMAT(stream,BF("x = %d y = %d\n") % x % y );
//  core::write_bf_stream(BF("x = %d y = %d\n") % x % y );
//  
//  if ( iColor != COLOR_DEPTH  && octTree->PfCharges.size() == 0 ) { //!octTree->PfCharges ) {
//    printf(( "charge coloring but no charges\n" ));
//    return;
//  }
//  fChargeFile = fopen(stream, "w"); //fChargeFile = fopen(sFileName, "w");
//	//if ( fChargeFile == NULL ) {
//  if (!fChargeFile) {
//    perror("file open error"); //perror(sFileName);
//    exit(1);
//  }
//  
//	/*
//	 *  Set up globals for this octree.
//	 */
//  
//  iColorMethod = iColor;
//  iMaxDepth = octTree->iMaxDepth;
//  dGridSize = octTree->dGridSize;
//  iNodeCount = 0;
//  
//	/*
//	 *  Print the nodes recursively.
//	 */
//  
//  octTree->OctNodePrintGrid( octTree->onHead );
//  
//  printf( "total nodes in octree %d\n", iNodeCount );
//  fclose( fChargeFile );
//  return;
//}


/*************************************************************************
 *************************************************************************/

/* global for speed */


//static void Octree_O::SplitIncludedNode( OctNode_sp PonNode )
void Octree_O::SplitIncludedNode( OctNode_sp PonNode)
{
  int		          i, j, k, nchild, ct, ct2, ccharge, cnode;
  gctools::Vec0<float>    PfTmpCharges;
  float                 PfCharge;
  int	boxct[8];
  
//	if ( PonNode->PonChildren != NULL )
//  for (i=0; i<8; i++){ 
//    if ( PonNode->PonChildren[i] && PonNode->PonChildren[i].notnilp())
//      printf("depth %d nodenum  %d status %d \n", PonNode->iDepth, PonNode->iNodeNum, PonNode->iStatus);
//      SIMPLE_ERROR(BF( "Programming error\n"));
//  }
	/*
	 *  Subdivide this node: set up children array
	 */
    PonMakeChildren(PonNode, PonNode->iDepth + 1, OCT_INCLUDED, this->iNodeCount );
//    printf("%s:%d:%s vPoint -> %d\n", __FILE__, __LINE__, __FUNCTION__,  PonNode->PonChildren[0]->iNodeNum);
  
	/*
	 *  copy the parent's charges to a temp array
	 */
  //GC_ALLOCATE(float, PfTmpCharges); // MALLOC( PfTmpCharges, float *, 
			//sizeof(float) * PiDensities[PonNode->iDepth]);
  //memcpy( PfTmpCharges, PonNode->PfCharges, 
  //        sizeof(float) * PiDensities[PonNode->iDepth]);
    

  for (i=0; i<PonNode->_PfCharges.size();i++){
    PfTmpCharges.push_back(PonNode->_PfCharges[i]);
  }

  for (i=0; i<this->_PfCharges.size(); i++){
    if (this->iNodeNumCharges[i] == PonNode->iNodeNum){
      cnode = i;
      break;
    }
  }
  
	/*
	 *  divide the parent's actual charge array space 
	 *	between the children (just assigning pointers)
	 */
  for (i=0; i<8; i++) {
    for (j=0; j<PonNode->_PfCharges.size(); j++){
//      PonNode->PonChildren[i]->_PfCharges.push_back(PonNode->_PfCharges[j+i * this->PiDensities[PonNode->iDepth+1]]);
      PonNode->PonChildren[i]->_PfCharges.push_back(this->_PfCharges[cnode+j+i * this->PiDensities[PonNode->iDepth+1]]);
    }
     }

	/*
	 *  Distribute charges to the children:
	 *	Reshuffle from the order in the current node
	 *	to the orders in the subboxes. I.e. while scanning
	 *	the copy of the parent's charge array in its order,
	 *	distribute the charges to the proper places in the 
	 *	childrens' arrays. 
	 *
	 *  For this to be relatively simple, the order of the children
	 *	in space (in the octree definition) mimics the order of 
	 *	evaluation of the parent's charges - see PonMakeChildren():
	 *
	 *		inner loop on Z axis (flagged by low order bit = '1')
	 *		next loop on Y axis (second bit = '2') 
	 *		outer loop on X axis (3rd bit = '4').
	 *
	 *	Use of the bits as described orders the subboxes properly.
	 */


	/*
	 *  initialize per-child point counts for indexing
	 *	childrens' point arrays
	 */
  for (i=0; i<8; i++)
    boxct[i] = 0;
  
	/*
	 *  nchild indicates child that current point of parent
	 *	will go to; 1st child has all bits off (000)
	 */
  nchild = 0;
  
	/*
	 *  total points in parent node
	 */
  ct = this->_iMaxDepth - PonNode->iDepth + 1;
  ct2 = ct / 2;
  
	/*
	 *  loop over points in parent node
	 */

  ccharge = 0;
  for (i=0; i<ct; i++) {
		/* loop over X axis */
    if (i==ct2)		/* 2nd half of X axis */
      nchild |= 4;	/* turn on X bit for 2nd slice*/
    nchild &= ~2;	/* turn off Y bit */
    for (j=0; j<ct; j++) {
			/* loop over Y axis */
      if (j==ct2)		/* 2nd half of Y axis */
        nchild |= 2;	/* Y bit on for 2nd slice */
      nchild &= ~1;	/* turn off Z bit */
      for (k=0; k<ct; k++) {
        PfCharge = PfTmpCharges[ccharge];
				/* loop over Z axis */
        if (k==ct2)		/* 2nd half of Z axis */
          nchild |= 1;	/* Z bit on for 2nd slice */

        PonNode->PonChildren[nchild]->_PfCharges[boxct[nchild]++] = PfCharge;

        ccharge++;
      }
    }
  }

	//FREE( PfTmpCharges );
}
//static int Octree_O::OctNodeDeleteSphere( OctNode_sp PonNode )
int Octree_O::OctNodeDeleteSphere( OctNode_sp PonNode,  double dDeleteRadius )
{
  int		i, iIncluded, iExcluded;
  double	d, dHalfEdge, dHalfDiagonal;
  OctNode_sp	PonChildren[8];
  Vector3	vCenter;
//
//  printf("nodestatus: %d iIncluded %d\n", PonNode->iStatus, iIncluded );
//  printf("num: %d iDepth %d\n", PonNode->iNodeNum, PonNode->iDepth );
 
	/*
	 *  If already excluded, quit. 
	 */
 
  if ( PonNode->iStatus == OCT_EXCLUDED ){
    return(OCT_EXCLUDED);
  }
	/*
	 *  If resolution reached, check corner only.
	 */
  if ( PonNode->iDepth == this->_iMaxDepth ) {
    d = dDistance( PonNode->vCorner, this->vNewPoint );
    if ( d < dDeleteRadius ) {
#ifdef OCTDEBUG
      printf(subtractions, ".dot %f %f %f\n", 
             PonNode->vCorner.getX(), 
             PonNode->vCorner.getY(), 
             PonNode->vCorner.getZ());
#endif
      PonNode->iStatus = OCT_EXCLUDED;
    }

    return(PonNode->iStatus);
  }
	/*
	 *  consider whole box from center
	 */
  dHalfEdge = this->PdHalfEdges[PonNode->iDepth];
  vCenter = PonNode->vCorner;
  vCenter.getX() += dHalfEdge;
  vCenter.getY() += dHalfEdge;
  vCenter.getZ() += dHalfEdge;
  d = dDistance( vCenter, this->vNewPoint );
  dHalfDiagonal = this->PdHalfDiagonals[PonNode->iDepth];

	/*
	 *  if box completely outside sphere, leave it alone.
	 */

  if ( d - dHalfDiagonal > dDeleteRadius ){
    return(OCT_INCLUDED);
  }

	/*
	 *  if box completely inside sphere, clean up its contents.
	 */
  if ( d + dHalfDiagonal < dDeleteRadius ) {
#ifdef SUB_DBG
    XYZ	point;
    int     j, k, ct = this->_iMaxDepth - oct->depth + 1;
    
    point.X = corner->X;
    for (i=0; i<ct; i++) {
      point.Y = corner->Y;
      for (j=0; j<ct; j++) {
        point.Z = corner->Z;
        for (k=0; k<ct; k++) {
          printf(subtractions, ".dot %f %f %f\n",
                 point.X, point.Y, point.Z);
          point.Z += this->dGridSize;
        }
        point.Y += this->dGridSize;
      }
      point.X += this->dGridSize;
    }
#endif
    //DestroyOctant( PonNode, OCT_EXCLUDED );
    PonNode->iStatus = OCT_EXCLUDED;

    return(OCT_EXCLUDED);
  }

	/*
	 *  recurse on sub-boxes
	 */

  if ( PonNode->iStatus == OCT_INCLUDED )
    SplitIncludedNode( PonNode);
  
	/*
	 *  Evaluate children.
	 */
//   printf("%s:%d:%s Entered\n", __FILE__, __LINE__, __FUNCTION__ );
 
  iIncluded = 0;
  iExcluded = 0;
  for (i=0; i<8; i++) {
      PonChildren[i] = PonNode->PonChildren[i];

      switch ( OctNodeDeleteSphere(PonChildren[i], dDeleteRadius ) ) {
      case OCT_INCLUDED:
          iIncluded++;	/* for debug tracking */
          break;
      case OCT_EXCLUDED:
          iExcluded++;
          break;
      case OCT_PARTIAL:
          break;
      default:
          SIMPLE_ERROR(BF("bad type\n"));
      }
//        perror("bad type\n");
//        exit(1);
  }

  if ( iExcluded == 8 ) {
    
		/*
		 *  Whole box is within vdw of multiple atoms.
		 */
    //DestroyOctant( PonNode, OCT_EXCLUDED );
    PonNode->iStatus = OCT_EXCLUDED;
//     printf("%s:%d:%s vPoint -> %d, %d, %d\n", __FILE__, __LINE__, __FUNCTION__, PonNode->iStatus, PonNode->iNodeNum, PonNode->iDepth);
//      printf("%s:%d:%s vPoint -> %d, %d, %d, %d\n", __FILE__, __LINE__, __FUNCTION__, PonNode->iStatus, PonNode->iNodeNum, PonNode->iDepth, i);

    return(OCT_EXCLUDED);
  }

	/*
	 * (don't fold into 1 node if iIncluded == 8 to avoid complications 
	 *	in regrouping charges from subnodes - so even if technically 
	 *	OCT_INCLUDED, it is evaluated as an OCT_PARTIAL)
	 */
//  printf("%s:%d:%s vPoint -> %d, %d, %d\n", __FILE__, __LINE__, __FUNCTION__, PonNode->iStatus, PonNode->iNodeNum, PonNode->iDepth);

  PonNode->iStatus = OCT_PARTIAL;
  return(OCT_PARTIAL);
}
CL_DEFMETHOD void Octree_O::OctTreeDeleteSphere( /*Octree_sp octTree,*/ Vector3 vPoint, double dRadius )
{
  
	/*
	 *  Set up globals for this octree.
	 */
  
//  PiDensities = octTree->PiDensities;
  this->vNewPoint = vPoint;
//  printf("%s:%d:%s vPoint -> %ld, %ld, %ld\n", __FILE__, __LINE__, __FUNCTION__, vPoint.getX(), vPoint.getY(), vPoint.getZ());
  double dDeleteRadius = dRadius;
//  dGridSize = octTree->dGridSize;
  
	/*
	 *  Delete recursively.
	 */

  OctNodeDeleteSphere( this->onHead,  dDeleteRadius  );
  return;
}



/*************************************************************************
 *************************************************************************/

//static void Octree_O::OctNodeUpdateCharge( OctNode_sp *PonNode, int iParentAtoms, Atom_sp *PaParentAtoms )
CL_DEFMETHOD void Octree_O::OctNodeUpdateCharge( OctNode_sp PonNode, int iParentAtoms, gctools::Vec0<Atom_sp> PaParentAtoms,  int iDistanceCharge )
{
  Vector3         	vPoint;
  double          	d;
  int             	i, j, k, l, ct, ccharge;
  int              	iCompOk;
  float            	PfCharge;
  Atom_sp	                PaAtom;
  gctools::Vec0<Atom_sp>  PaAtoms;
  int             	iAtoms;

 

	/*
	 *  If partial, recurse.
	 */
  if ( PonNode->iStatus == OCT_PARTIAL ) {
    for (i=0; i<8; i++)
      this->OctNodeUpdateCharge( PonNode->PonChildren[i],
                                 PonNode->iAtoms, PonNode->PaAtomList, iDistanceCharge );
    return;
  }

  
	/*
	 *  If already excluded, quit.
	 */
  if ( PonNode->iStatus == OCT_EXCLUDED )
    return;
  
    
	/*
	 *  Included, do whole thing
	 */
  
//	if ( PonNode->PaAtomList != NULL ) {
  if ( PonNode->PaAtomList.size() > 0 ){
    for (i=0; i<PonNode->PaAtomList.size(); i++){      
      PaAtoms.push_back(PonNode->PaAtomList[i]);
    }
    iAtoms = PonNode->iAtoms;
  } else {
    for (i=0; i<PaParentAtoms.size(); i++){      
      PaAtoms.push_back(PaParentAtoms[i]);
    }
    iAtoms = iParentAtoms;
  }
  
  ct = this->_iMaxDepth - PonNode->iDepth + 1;
  ccharge = 0;
  vPoint.getX() = PonNode->vCorner.getX();
  for (i=0; i<ct; i++, vPoint.getX()+=this->dGridSize) {
    vPoint.getY() = PonNode->vCorner.getY();
    for (j=0; j<ct; j++, vPoint.getY()+=this->dGridSize) {
      vPoint.getZ() = PonNode->vCorner.getZ();
      for (k=0; k<ct; k++, vPoint.getZ()+=this->dGridSize) {
        PfCharge = PonNode->_PfCharges[ccharge];
        double 	dX, dY, dZ;

        
        dX = this->vNewPoint.getX() - vPoint.getX();
        dX = dX * dX;
        dY = this->vNewPoint.getY() - vPoint.getY();
        dY = dY * dY;
        dZ = this->vNewPoint.getZ() - vPoint.getZ();
        dZ = dZ * dZ;
        d = dX + dY + dZ;
        if ( iDistanceCharge )
          d = sqrt(d);

        PfCharge += this->fNewCharge / d;
			/*
			 *  Look at neighboring atoms to make sure
			 *	this point is valid for min/max check.
			 */
        iCompOk = 1;
        for (l=0; l<iAtoms; l++) {
        PaAtom = PaAtoms[l];
                          //dX = vAtomPosition( *PaAtom ).getX() - vPoint.getX();
          dX = PaAtom->getPosition().getX() - vPoint.getX();
          dX = dX * dX;
                         // dY = vAtomPosition( *PaAtom ).getY() - vPoint.getY();
          dY = PaAtom->getPosition().getY() - vPoint.getY();
          dY = dY * dY;
                         // dZ = vAtomPosition( *PaAtom ).getZ() - vPoint.getZ();
          dZ = PaAtom->getPosition().getZ() - vPoint.getZ();
          dZ = dZ * dZ;
          d = dX + dY + dZ;
                          //if ( d < AtomTemp(PaAtom) ) {
          if ( d < PaAtom->dAtomTemp ) {
            iCompOk = 0;
            break;
          }
        }

        if ( iCompOk ) {
          if (PfCharge > this->fMaxCharge ) {
            this->fMaxCharge = PfCharge;
            this->vMaxCharge = vPoint;
          } else if (PfCharge < this->fMinCharge ) {
            this->fMinCharge = PfCharge;
            this->vMinCharge = vPoint;
          }
        } else {
				/* HACK to ensure printgrid coloring ok */
          PfCharge = 0.0;
        }
//        printf("%s:%d:%s vPoint -> %d, %d, %d, %d, %d, %d, %f\n", __FILE__, __LINE__, __FUNCTION__, PonNode->iStatus, PonNode->iNodeNum, PonNode->iDepth, i, j, k, PfCharge);

        PonNode->_PfCharges[ccharge] = PfCharge;
        ccharge++;
      }
    }
  }
  
  for (l=0; l< PonNode->_PfCharges.size(); l++) {
    this->_PfCharges.push_back(PonNode->_PfCharges[l]);
    this->iNodeNumCharges.push_back(PonNode->iNodeNum);
  }

}
//CL_DEFMETHOD void Octree_O::OctTreeUpdateCharge( Octree_sp octTree, Vector3 vNewPoint, float fCharge, double dCutDist,
//                          Vector3 vMax, Vector3 vMin )
CL_DEFMETHOD core::T_mv Octree_O::OctTreeUpdateCharge( /*Octree_sp octTree,*/ Vector3 vNewPoint, float fCharge, double dCutDist)
{
  Atom_sp	PaAtom;
  int	i, iDistanceCharge;
  int iChargeAtoms = this->vaAtoms.size();

  if ( this->type != Shell ) {
    SIMPLE_ERROR(BF("UpdateCharge: wrong tree type\n"));
  }
  if ( this->_PfCharges.size()==0 ) {//( !octTree->PfCharges ) {
    SIMPLE_ERROR(BF("UpdateCharge: charges not initted\n"));
  }
  
	/*
	 *  Set up globals for octree.
	 */
  this->vNewPoint = vNewPoint;
  this->fNewCharge = fCharge;
  this->fMaxCharge = -FLT_MAX;
  this->fMinCharge = FLT_MAX;
//  dGridSize = octTree->dGridSize;
  if ( this->iDielectric == DIEL_R2 )
    iDistanceCharge = 0;
  else
    iDistanceCharge = 1;
  
	/*
	 *  Add new radius to atoms.
	 */
	//PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
  for (i=0; i<iChargeAtoms; i++){ // , PaAtom++) {
    this->vaAtoms[i]->dAtomTemp += dCutDist; // AtomTempDoubleIncrement( *PaAtom, dCutDist );
    this->vaAtoms[i]->dAtomTemp *= this->vaAtoms[i]->dAtomTemp; //AtomTempDoubleSquare( *PaAtom );
  }

  this->_PfCharges.clear();
  this->iNodeNumCharges.clear();

	/*
	 *  Descend octree.
	 */
  this->OctNodeUpdateCharge( this->onHead, iChargeAtoms, 
                                this->vaAtoms, iDistanceCharge) ; //PVAI( octTree->vaAtoms, ATOM, 0 ) );
//  vMin = vMinCharge;
//  vMax = vMaxCharge;
	/*
	 *  Restore atoms' radii.
	 */
	//PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
  for (i=0; i<iChargeAtoms; i++ ){ //, PaAtom++) {
    this->vaAtoms[i]->dAtomTemp = sqrt(this->vaAtoms[i]->dAtomTemp); // AtomTempDoubleSquareRoot( *PaAtom );
    this->vaAtoms[i]->dAtomTemp += -dCutDist; // AtomTempDoubleIncrement( *PaAtom, -dCutDist );
  }
    return Values(geom::OVector3_O::create(this->vMinCharge),geom::OVector3_O::create(this->vMaxCharge));
//  return;
}


/*************************************************************************
 *************************************************************************/
/*
 *  Solvent check: see if a solvent overlaps a new ion. 
 *	Descend tree nodes which contain the new point inside their
 *	box-including sphere till a final node is reached, then
 *	check the atom list for that node.
 *
 *  NOTE - this code assumes that atom lists are updated as solvent
 *	is deleted. In fact, this is not done, so this code is left
 *	for adaptation in finding vdw pairs.
 */

//static int Octree_O::OctNodeCheckSolvent( OctNode_sp PonNode )
int Octree_O::OctNodeCheckSolvent( OctNode_sp PonNode )
{
  int	i;
  Atom_sp	PaAtom;
  Vector3	vCenter;
  double	d, dHalfEdge, dHalfDiagonal;
  bool          isPonChildren;
  
  if ( PonNode->iStatus != OCT_INCLUDED )
    return(0);
  
	/*
	 *  See if newpoint is in minimally box-enclosing sphere
	 */
  
  dHalfEdge = this->PdHalfEdges[PonNode->iDepth];
  vCenter = PonNode->vCorner;
  vCenter.getX() += dHalfEdge;
  vCenter.getY() += dHalfEdge;
  vCenter.getZ() += dHalfEdge;
  d = dDistance( vCenter, this->vNewPoint );
  dHalfDiagonal = this->PdHalfDiagonals[PonNode->iDepth];
  
  if ( d > dHalfDiagonal ) {
    return(0);
  }
/*
VP0(("lev %d within %f %f %f  %f %f %f\n", 
PonNode->iDepth,
PonNode->vCorner.dX,
PonNode->vCorner.dY,
PonNode->vCorner.dZ,
PonNode->vCorner.dX + 2 * PdHalfEdges[PonNode->iDepth],
PonNode->vCorner.dY + 2 * PdHalfEdges[PonNode->iDepth],
PonNode->vCorner.dZ + 2 * PdHalfEdges[PonNode->iDepth]));
*/
	/*
	 *  If at last depth w/ atomlist, check them.
	 */
//	if ( PonNode->PonChildren == NULL ) {
  isPonChildren = false;
  for (i=0; i<8; i++){
    if ( PonNode->PonChildren[i].notnilp()){
      isPonChildren = true;
    }
  }
  
  if(!isPonChildren)
  {
/*
VP0(("final? \n"));
*/
    for (i=0; i<PonNode->iAtoms; i++) {
      PaAtom = PonNode->PaAtomList[i];
      d = dDistanceSq( this->vNewPoint,
                       PaAtom->getPosition() );
/*
VP0(("d= %f \n", d));
*/
      if ( d < dClosestDistance ) {
        this->aClosestAtom = PaAtom;
        this->dClosestDistance = d;
      }
    }
    return(1);
  }
  
	/*
	 *  partial: subdivide
	 */
  for (i=0; i<8; i++) 
    if( OctNodeCheckSolvent( PonNode->PonChildren[i] ) )
      return(1);
  
  return(0);
}



Residue_sp Octree_O::rOctTreeCheckSolvent( /*Octree_sp octTree,*/ Vector3 vPoint )
{
  //Atom_sp	PaAtom;
  if ( this->type != InteriorSolvent ) {
    SIMPLE_ERROR(BF("CheckSolvent: wrong octree type\n"));
//    printf( "CheckSolvent: wrong octree type\n" );
//    exit(1);
  }
	/*
	 *  Set up globals for octree.
	 */
  this->vNewPoint = vPoint;
//  dGridSize = octTree->dGridSize;
  
	/*
	 *  Initialize closest atom for comparison.
	 */
  
	//PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 ); 
  //PaAtom = octTree->vaAtoms[0]; 
  this->aClosestAtom = this->vaAtoms[0]; //*PaAtom;
  this->dClosestDistance = this->dDistanceSq( this->vNewPoint,
                                  this->aClosestAtom->getPosition());
  
	/*
	 *  Descend octree.
	 */
  if (!this->OctNodeCheckSolvent( this->onHead ) ) {
    printf("Completely out of solvent bounding area\n");
//		return(NULL);
    return(_Nil<core::T_O>());
  }
  
	/*
	 *  Check if closest atom overlaps at all.
	 */
  if ( this->dClosestDistance < 9.0 ) {	/* HACK - approx ion+wat */
    vPoint = this->vNewPoint;
    return gc::As<Residue_sp>(this->aClosestAtom->containedBy());
    //return gc::As_unsafe<Residue_sp>(aClosestAtom->containedBy());
    // return( (Residue_sp) cContainerWithin( aClosestAtom ) );
  }
  printf("No overlap w/ solvent\n");
//	return(NULL);
  return(_Nil<core::T_O>());
}

};
