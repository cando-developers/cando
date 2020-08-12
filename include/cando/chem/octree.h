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
#include <clasp/core/array.h>
#include <cando/chem/atom.h>
#include <cando/geom/vector3.h>
#include <cando/chem/chemPackage.h>
#include <cando/chem/ffNonbondDb.h>

namespace chem {
FORWARD(BoundingBox);
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



FORWARD(AddIonOctree);
class AddIonOctree_O : public core::CxxObject_O {
  LISP_CLASS(chem,ChemPkg,AddIonOctree_O,"AddIonOctree",core::CxxObject_O);
public:
    bool fieldsp() const { return true;};
    void fields(core::Record_sp record);
  typedef enum { Shell, InteriorSolute, InteriorSolvent } AddIonOctreeEnum;
public:
  AddIonOctreeEnum                    type;
  BoundingBox_sp                _BoundingBox;
  double                        _BoundingBox_x_rsize;
  double                        _BoundingBox_y_rsize;
  double                        _BoundingBox_z_rsize;
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
  core::HashTableEq_sp          atomsToResidues;
  core::HashTableEq_sp          residuesToMolecules;
  OctNode_sp                    onHead;

public:
  static AddIonOctree_sp make(Aggregate_sp uUnit, AddIonOctreeEnum iType, double dGridSpace, 
                         double dAddExtent, double dShellExtent, FFNonbondDb_sp nonbondDb, int iIncludeSolvent, bool bVerbose);
public:
  void setBoundingBox(BoundingBox_sp bounding_box);
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
  int dFinalCheck( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp>& PaAtomList, double dShellRadius  );
  void DestroyOctant( OctNode_sp PonNode, int iStatus );
  int iBuildShellOctant( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp>& PaAtomList,  double dShellRadius );
  int iBuildInteriorOctant( OctNode_sp PonNode, int iAtoms, gctools::Vec0<Atom_sp>& PaAtomList );
  void OctreeDestroy();
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
  void OctNodeUpdateCharge( OctNode_sp PonNode, int iParentAtoms, gctools::Vec0<Atom_sp>& PaParentAtoms, int iDistanceCharge );
  int OctNodeCheckSolvent( OctNode_sp PonNode );
  void OctreeCreate(Aggregate_sp uUnit, AddIonOctreeEnum iType, double dGridSpace, 
                     double dAddExtent, double dShellExtent, FFNonbondDb_sp nonbondDb, int iIncludeSolvent, bool bVerbose);
  core::T_mv OctreeInitCharges(/*AddIonOctree_sp octTree,*/ int iAtomOption, int iDielectric, 
                           double dCutDist);	
  void OctreeDescribe();	/* ( AddIonOctree_sp ) */
  void OctreeDeleteSphere(/*AddIonOctree_sp octTree, */ Vector3 vPoint, double dRadius);
  core::T_mv OctreeUpdateCharge(/*AddIonOctree_sp octTree,*/ Vector3 vNewPoint, 
                            float fCharge, double dCutDist);
  core::T_sp rOctreeCheckSolvent(/*AddIonOctree_sp octTree,*/ Vector3 vPoint);

  AddIonOctree_O() : _BoundingBox(_Unbound<BoundingBox_O>()) {};
//  void OctreePrintGrid(AddIonOctree_sp octTree, core::T_sp stream, int iColor);

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

DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::AddIonOctree_O::AddIonOctreeEnum, chem::_sym_STARoctree_typeSTAR);



/**

Based on code from: https://github.com/brandonpelfrey/SimpleOctree

MIT License

Copyright (c) 2017 Brandon Pelfrey

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*!
 * A simple octree implementation based on Brandon Pelfrey's code.
 */

namespace chem {
FORWARD(GenericOctree);
};
template <>
struct gctools::GCInfo<chem::GenericOctree_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};


namespace chem {
FORWARD(GenericOctree);
class GenericOctree_O : public core::CxxObject_O {
  LISP_CLASS(chem,ChemPkg,GenericOctree_O,"GenericOctree",core::CxxObject_O);
		// Physical position/size. This implicitly defines the bounding 
		// box of this node
  size_t  _depth;
  Vector3 _origin;         //! The physical center of this node
  Vector3 _halfDimension;  //! Half the width/height/depth of this node
  float   _containing_radius;
  const size_t Entries = 2;
		// The tree has up to eight children and can additionally store
		// a point, though in many applications only, the leaves will store data.
  GenericOctree_sp _children[8]; //! Pointers to child octants
  core::T_sp   _data;   //! Data points to be stored at a node
  Vector3  _position; //! Positions of data stored at a node

		/*
				Children follow a predictable pattern to make accesses simple.
				Here, - means less than 'origin' in that dimension, + means greater than.
				child:	0 1 2 3 4 5 6 7
				x:      - - - - + + + +
				y:      - - + + - - + +
				z:      - + - + - + - +
		 */
public:
  bool fieldsp() const { return true;};
  void fields(core::Record_sp record);
public:
  CL_DOCSTRING(R"doc(Make a generic-octree centered on ORIGIN with HALF-DIMENSIONs)doc");
  CL_LISPIFY_NAME(make-generic-octree);
  CL_DEF_CLASS_METHOD
  static GenericOctree_sp make(const Vector3& origin, const Vector3& halfDimension);
public:
  GenericOctree_O(const Vector3& origin, const Vector3& halfDimension) 
    : _depth(0), _origin(origin), _halfDimension(halfDimension), _data(_Unbound<core::T_O>()) {
				// Initially, there are no children
    for(int i=0; i<8; ++i) 
      this->_children[i] = _Unbound<GenericOctree_O>();
    this->_containing_radius = this->_halfDimension.getX()*this->_halfDimension.getX()
      + this->_halfDimension.getY()*this->_halfDimension.getY()
      + this->_halfDimension.getZ()*this->_halfDimension.getZ();
  }


  size_t depth() const;
  Vector3 origin() const;
  Vector3 halfDimension() const;
  float   containing_radius() const;
  GenericOctree_mv children() const;
  core::T_sp data() const;
  bool dataBoundP() const;
  Vector3 position() const;
  
		// Determine which octant of the tree would contain 'point'
  int getOctantContainingPoint(const Vector3& point) const {
    int oct = 0;
    if(point.getX() >= _origin.getX()) oct |= 4;
    if(point.getY() >= _origin.getY()) oct |= 2;
    if(point.getZ() >= _origin.getZ()) oct |= 1;
    return oct;
  }

  bool isLeafNode() const;

  void insert_impl(const Vector3& point, core::T_sp data) {
			// If this node doesn't have a data point yet assigned 
			// and it is a leaf, then we're done!
    if (this->_depth>20) {
      SIMPLE_ERROR(BF("The octree is too deep - 20 levels"));
    }
    if(isLeafNode()) {
      if (this->_data.unboundp()) {
        this->_data = data;
        this->_position = point;
        return;
      } else {
					// We're at a leaf, but there's already something here
					// We will split this node so that it has 8 child octants
					// and then insert the old data that was here, along with 
					// this new data point

					// Save this data point that was here for a later re-insert
        core::T_sp oldData = this->_data;
        this->_data = _Unbound<core::T_O>();
        Vector3 oldPosition = this->_position;
        double dx = fabs(point.getX()-oldPosition.getX());
        double dy = fabs(point.getY()-oldPosition.getY());
        double dz = fabs(point.getZ()-oldPosition.getZ());
        double dist_squared = dx*dx+dy*dy+dz*dz;
        if (dist_squared<0.1) {
          SIMPLE_ERROR(BF("oldPosition %s and point %s are too close to add to octree") % oldPosition.asString() % point.asString());
        }

					// Split the current node and create new empty trees for each
					// child octant.
        for(int i=0; i<8; ++i) {
						// Compute new bounding box for this child
          Vector3 new_Origin = this->_origin;
          new_Origin.getX() += _halfDimension.getX() * (i&4 ? .5f : -.5f);
          new_Origin.getY() += _halfDimension.getY() * (i&2 ? .5f : -.5f);
          new_Origin.getZ() += _halfDimension.getZ() * (i&1 ? .5f : -.5f);
          GC_ALLOCATE_VARIADIC(GenericOctree_O,child,new_Origin, this->_halfDimension*.5f);
          child->_depth = this->_depth+1;
          this->_children[i] = child;
        }

					// Re-insert the old point, and insert this new point
					// (We wouldn't need to insert from the root, because we already
					// know it's guaranteed to be in this section of the tree)
        this->_children[getOctantContainingPoint(oldPosition)]->insert_impl(oldPosition,oldData);
        this->_children[getOctantContainingPoint(point)]->insert_impl(point,data);
      }
    } else {
				// We are at an interior node. Insert recursively into the 
				// appropriate child octant
      int octant = getOctantContainingPoint(point);
      this->_children[octant]->insert_impl(point,data);
    }
  }


  CL_LISPIFY_NAME(generic-octree-insert);
  CL_DEFMETHOD
  void insert(const Vector3& point, core::T_sp data) {
			// If this node doesn't have a data point yet assigned 
			// and it is a leaf, then we're done!
    if (this->_depth==0 &&
        ((point.getX()<(this->_origin.getX()-this->_halfDimension.getX())) ||
         (point.getY()<(this->_origin.getY()-this->_halfDimension.getY())) ||
         (point.getZ()<(this->_origin.getZ()-this->_halfDimension.getZ())) ||
         (point.getX()>(this->_origin.getX()+this->_halfDimension.getX())) ||
         (point.getY()>(this->_origin.getY()+this->_halfDimension.getY())) ||
         (point.getZ()>(this->_origin.getZ()+this->_halfDimension.getZ())))) {
      SIMPLE_ERROR(BF("Trying to insert out of bounds point %lf,%lf,%lf into octree min %lf,%lf,%lf  max %lf,%lf,%lf") % point.getX() % point.getY() % point.getZ()
                   % (this->_origin.getX()-this->_halfDimension.getX())
                   % (this->_origin.getY()-this->_halfDimension.getY())
                   % (this->_origin.getZ()-this->_halfDimension.getZ())
                   % (this->_origin.getX()+this->_halfDimension.getX())
                   % (this->_origin.getY()+this->_halfDimension.getY())
                   % (this->_origin.getZ()+this->_halfDimension.getZ()));
    }
    this->insert_impl(point,data);
  }

		// This is a really simple routine for querying the tree for points
		// within a distance of the querypoint inside of a periodic box
		// All results are pushed into 'results'
  void getPointsWithinCutoff(double cutoff_squared, double cutoff, core::T_sp query_data, const Vector3& querypoint, double x_size, double y_size, double z_size, double x_rsize, double y_rsize, double z_rsize, const Matrix& from_target_transform, core::ComplexVector_sp results );
  void getPointsWithinCutoffNoTransform(double cutoff_squared, double cutoff, core::T_sp query_data, const Vector3& querypoint, double x_size, double y_size, double z_size, double x_rsize, double y_rsize, double z_rsize, core::ComplexVector_sp results );
  void walkGenericOctree(core::Function_sp callback);
  bool octreeTransformedIntersectsInBoundingBox(const Matrix& mytransform, GenericOctree_sp other, const Matrix& otherTransform, BoundingBox_sp boundingBox);
};

};


#endif
