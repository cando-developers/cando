#define	DEBUG_LEVEL_FULL


//
//	superpose
//
//	Superpose two sets of points,
//	return the transformation that transforms moveable
//	to the fixed ones.
//




#include <stdio.h>
#include <math.h>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/geom/matrix.h>
#include <cando/geom/omatrix.h>
#include <cando/chem/symbolTable.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/loop.h>
#include <cando/chem/superposeEngine.h>
#include <clasp/core/intArray.h>
#include <clasp/core/wrappers.h>


namespace chem
{

    SYMBOL_EXPORT_SC_(ChemPkg,superpose);


    EXPOSE_CLASS_AND_GLOBALS(chem,SuperposeEngine_O);
    EXPOSE_CLASS(chem,SuperposeSelectedAtoms_O);

    void SuperposeEngine_O::exposeCando(core::Lisp_sp e)
    {
	core::class_<SuperposeEngine_O>()
	    .def("setFixedPoints",&SuperposeEngine_O::setFixedPoints)
	    .def("setFixedAllPoints",&SuperposeEngine_O::setFixedAllPoints)
	    .def("setMoveablePoints",&SuperposeEngine_O::setMoveablePoints)
	    .def("setMoveableAllPoints",&SuperposeEngine_O::setMoveableAllPoints)
	    .def("getNumberOfFixedPoints",&SuperposeEngine_O::getNumberOfFixedPoints)
	    .def("getNumberOfMoveablePoints",&SuperposeEngine_O::getNumberOfMoveablePoints)
	    .def("superpose",&SuperposeEngine_O::superpose)
	    .def("rootMeanSquareDifference",&SuperposeEngine_O::rootMeanSquareDifference)
	    .def("debugString",&SuperposeEngine_O::debugString)
	    ;
    }

    void SuperposeEngine_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,SuperposeEngine,"","",_lisp)
	    .def("setFixedPoints",&SuperposeEngine_O::setFixedPoints)
	    .def("setFixedAllPoints",&SuperposeEngine_O::setFixedAllPoints)
	    .def("setMoveablePoints",&SuperposeEngine_O::setMoveablePoints)
	    .def("setMoveableAllPoints",&SuperposeEngine_O::setMoveableAllPoints)
	    .def("getNumberOfFixedPoints",&SuperposeEngine_O::getNumberOfFixedPoints)
	    .def("getNumberOfMoveablePoints",&SuperposeEngine_O::getNumberOfMoveablePoints)
	    .def("superpose",&SuperposeEngine_O::superpose)
	    .def("rootMeanSquareDifference",&SuperposeEngine_O::rootMeanSquareDifference)
	    .def("debugString",&SuperposeEngine_O::debugString)
	    ;
#endif
    }



    void SuperposeEngine_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();
    }



#ifdef XML_ARCHIVE
    void	SuperposeEngine_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("fixedIndices",this->_FixedIndices);
    node->attribute("fixedCoordinates",this->_FixedCoordinates);
    node->attribute("moveableIndices",this->_MoveableIndices);
    node->attribute("moveableCoordinates",this->_MoveableCoordinates);
}
#endif



    string SuperposeEngine_O::debugString()
    {
	stringstream ss;
	for (uint i=0; i<this->_FixedIndices->size(); i++ )
	{
	    int idx = this->_FixedIndices->get(i);
	    ss << "Fixed#" << i << " " << this->_FixedCoordinates->getElement(idx).asString() << std::endl;
	}
	for (uint j=0; j<this->_MoveableIndices->size(); j++ )
	{
	    int idx = this->_MoveableIndices->get(j);
	    ss << "Moveable#" << j << " " << this->_MoveableCoordinates->getElement(idx).asString() << std::endl;
	}
	return ss.str();
    }









//
//	geometricCenterOfPoints
//
    Vector3	geometricCenterOfPointsIndirect( core::IntArray_sp indices, geom::CoordinateArray_sp pnts)
{ _G();
    Vector3				pos;
    core::IntArray_O::iterator		it;

    ASSERTP(indices->size()<=pnts->size(),"There must be at least as many indices as coordinates");
    pos.set(0.0,0.0,0.0);
    if ( indices->size() == 0 ) return pos;

    for ( it=indices->begin();
	  it!=indices->end();
	  it++ ) {
	pos = pos + pnts->getElement(*it);
    }
    pos = pos.multiplyByScalar(1.0/((double)(indices->size())));
    return pos;
}

void SuperposeEngine_O::initialize()
{_G();
    this->Base::initialize();
    this->_FixedIndices = core::IntArray_O::create();
    this->_MoveableIndices = core::IntArray_O::create();
    this->_FixedCoordinates = geom::CoordinateArray_O::create();
    this->_MoveableCoordinates = geom::CoordinateArray_O::create();
}



void	SuperposeEngine_O::eraseMoveablePoints()
{
    this->_MoveableIndices->clear();
    this->_MoveableCoordinates->clear();
}

void	SuperposeEngine_O::eraseFixedPoints()
{
    this->_FixedIndices->clear();
    this->_FixedCoordinates->clear();
}

    void SuperposeEngine_O::appendMoveablePoint(const Vector3& pos)
{_OF();
    ASSERT(this->_MoveableIndices.notnilp());
    ASSERT(this->_MoveableCoordinates.notnilp());
    LOG(BF("appendMoveablePoint vector: %s") % pos.asString() );
    int idx = this->_MoveableCoordinates->size();
    this->_MoveableIndices->append(idx);
    this->_MoveableCoordinates->push_back(pos);
    LOG(BF("There are now %d MoveableIndices") % this->_MoveableIndices->size());
    LOG(BF("There are now %d MoveableCoordinates") % this->_MoveableCoordinates->size());
}

    void SuperposeEngine_O::appendFixedPoint(const Vector3& pos)
{_OF();
    ASSERT(this->_FixedIndices.notnilp());
    ASSERT(this->_FixedCoordinates.notnilp());
    LOG(BF("appendFixedPoint vector: %s") % pos.asString() );
    int idx = this->_FixedCoordinates->size();
    this->_FixedIndices->append(idx);
    this->_FixedCoordinates->push_back(pos);
    LOG(BF("There are now %d FixedIndices") % this->_FixedIndices->size());
    LOG(BF("There are now %d FixedCoordinates") % this->_FixedCoordinates->size());
}



int	SuperposeEngine_O::getNumberOfFixedPoints()
{_G();
    return this->_FixedIndices->size();
}


int	SuperposeEngine_O::getNumberOfMoveablePoints()
{_G();
    return this->_MoveableIndices->size();
}


void	SuperposeEngine_O::doSuperpose()
{_G();
    VectorVector3s			Sj;
    VectorVector3s			Si;
    VectorVector3s::iterator	itS,itSi,itSj;
    core::IntArray_O::iterator		iaV;
    vector<VectorVector3s::iterator>::iterator	ititV;
    Vector3				fixedCenter,vTemp;
    Vector3				moveableCenter;
    Matrix				mat,M, MT, trM, P,evecs,em,trans,rot;
    Vector4				v4,evals,largestEv,quaternion;
    int				iLargestEv;
    double				X0, X1, X2, x0, x1, x2;
    double				l,m,n,s,ll,mm,nn,ss;
    int				fixedIndicesSize, moveableIndicesSize;
    fixedIndicesSize = this->_FixedIndices->size();
    moveableIndicesSize = this->_MoveableIndices->size();
    LOG(BF("Moveable indices fixed(%d) moveable(%d)") % fixedIndicesSize % moveableIndicesSize );
    ASSERTP(fixedIndicesSize==moveableIndicesSize,"num. fixed points must equal num. of moveable points");
    ASSERTF(fixedIndicesSize>=3, BF("You must have at least 3 points to superpose and you only have: %d")% fixedIndicesSize );
    Sj.resize(this->_FixedIndices->size());
    if ( this->_MoveableIndices->size() < 3 )
    {
	SIMPLE_ERROR(BF("Number of MoveableIndices must be greater than 3"));
    }
    Si.resize(this->_MoveableIndices->size());
    LOG(BF("this->_FixedIndices->size()=%d") % this->_FixedIndices->size()  );
    LOG(BF("this->_MoveableIndices->size()=%d") % this->_MoveableIndices->size()  );

    {_BLOCK_TRACE("Calculating geometric center of fixed points");
	fixedCenter = geometricCenterOfPointsIndirect(this->_FixedIndices,this->_FixedCoordinates);
        LOG(BF( "Translating fixed to fixed geometric center: %s")% 
	    fixedCenter.asString().c_str() );
	for ( iaV=this->_FixedIndices->begin(),itS=Sj.begin();
	      iaV!=this->_FixedIndices->end();
	      iaV++,itS++ ) {
	    LOG(BF("index=%d   *iaV=%d    pnt=%s") % (iaV-this->_FixedIndices->begin()) 
		% *iaV % this->_FixedCoordinates->getElement(*iaV).asString().c_str() ); 
	    *itS = this->_FixedCoordinates->getElement(*iaV) - fixedCenter;
	}
    }
    { _BLOCK_TRACE("Calculating geometricCenterOfPoints of moveable points");
	moveableCenter.set(0.0,0.0,0.0);
	for ( iaV=this->_MoveableIndices->begin(); iaV !=this->_MoveableIndices->end(); iaV++ ) {
	    moveableCenter = this->_MoveableCoordinates->getElement(*iaV) + moveableCenter;
	}
	moveableCenter = moveableCenter.multiplyByScalar(1.0/(double)(this->_MoveableIndices->size()));
	LOG(BF("Translating moveable to moveable geometric center: %s") % moveableCenter.asString().c_str()  );
	for ( iaV=this->_MoveableIndices->begin(),itS=Si.begin();
	      itS!=Si.end();
	      iaV++,itS++ ) {
	    *itS = this->_MoveableCoordinates->getElement(*iaV) - moveableCenter;
	    LOG(BF("index=%d   *iaV=%d    pnt=%s")
		% (iaV-this->_MoveableIndices->begin())
		% *iaV
		% this->_FixedCoordinates->getElement(*iaV).asString().c_str() );
	}
    }
    LOG(BF( "fixedCenter = %s")% fixedCenter.asString() );
    LOG(BF( "moveableCenter = %s")% moveableCenter.asString() );
    itSj = Sj.begin();
    itSi = Si.begin();
    M.setValue(0.0);
    for ( uint c = 0; c< Sj.size(); c++ ) {
	X0 = itSj->getX();
	X1 = itSj->getY();
	X2 = itSj->getZ();
	x0 = itSi->getX();
	x1 = itSi->getY();
	x2 = itSi->getZ();
        LOG(BF("Iterating c=%d") % c  );
        LOG(BF("Iterating x0=%lf") % x0  );
        LOG(BF("Iterating x1=%lf") % x1  );
        LOG(BF("Iterating x2=%lf") % x2  );
        LOG(BF("Iterating X0=%lf") % X0  );
        LOG(BF("Iterating X1=%lf") % X1  );
        LOG(BF("Iterating X2=%lf") % X2  );
	M.atRowColPut( 0, 0, M.atRowCol( 0, 0 )+x0*X0 );
	M.atRowColPut( 0, 1, M.atRowCol( 0, 1 )+x0*X1 );
	M.atRowColPut( 0, 2, M.atRowCol( 0, 2 )+x0*X2 );
	M.atRowColPut( 1, 0, M.atRowCol( 1, 0 )+x1*X0 );
	M.atRowColPut( 1, 1, M.atRowCol( 1, 1 )+x1*X1 );
	M.atRowColPut( 1, 2, M.atRowCol( 1, 2 )+x1*X2 );
	M.atRowColPut( 2, 0, M.atRowCol( 2, 0 )+x2*X0 );
	M.atRowColPut( 2, 1, M.atRowCol( 2, 1 )+x2*X1 );
	M.atRowColPut( 2, 2, M.atRowCol( 2, 2 )+x2*X2 );
	itSj++;
	itSi++;
    }
    LOG(BF( "M= %s")% M.asString() );
    MT = M.transpose();
    v4.set( 0.0, 0.0, 0.0, M.trace() );
    trM.setFromQuaternion(v4);
    P = M + (MT - (trM*2.0));
    P.atRowColPut(0,3,M.atRowCol(1,2)-M.atRowCol(2,1));
    P.atRowColPut(3,0,P.atRowCol(0,3));
    P.atRowColPut(1,3,M.atRowCol(2,0)-M.atRowCol(0,2));
    P.atRowColPut(3,1,P.atRowCol(1,3));
    P.atRowColPut(2,3,M.atRowCol(0,1)-M.atRowCol(1,0));
    P.atRowColPut(3,2,P.atRowCol(2,3));
    P.atRowColPut(3,3,0.0);
    LOG(BF("P=%s")% P.asString() );
    P.eigenSystem( evals, evecs );
    iLargestEv = evals.indexOfLargestElement();
    largestEv = evecs.getCol(iLargestEv);
    LOG(BF( "largestEv=%s")% largestEv.asString() );
    em.setFromQuaternion(largestEv);

    LOG(BF( "eigenMatrix= %s")% em.asString() );
    quaternion = em.getCol(3);
    LOG(BF( "quaternion = %s")% quaternion.asString() );
    //
    // Convert the quaternion into a rotation matrix
    rot.setValue(0.0);
    l = quaternion.getW();
    m = quaternion.getX();
    n = quaternion.getY();
    s = quaternion.getZ();
    ll = l*l;
    mm = m*m;
    nn = n*n;
    ss = s*s;
    rot.atRowColPut( 0,0, ll - mm - nn + ss );
    rot.atRowColPut( 0,1,  2*(l*m - n*s));
    rot.atRowColPut( 0,2,  2*(l*n + m*s));
    rot.atRowColPut( 1,0,  2*(l*m + n*s));
    rot.atRowColPut( 1,1,  -ll + mm - nn + ss);
    rot.atRowColPut( 1,2,  2*(m*n - l*s));
    rot.atRowColPut( 2,0,  2*(l*n - m*s));
    rot.atRowColPut( 2,1,  2*(m*n + l*s));
    rot.atRowColPut( 2,2,  -ll - mm + nn + ss);
    rot.atRowColPut( 3,3,  1.0 );
    LOG(BF( "rot=%s")% rot.asString());
    vTemp = moveableCenter*-1.0;
    trans.translate(&vTemp);
    LOG(BF( "moveableTrans=%s")% trans.asString() );
//    mat.setFromQuaternion(quaternion);
    mat = rot*trans;
    LOG(BF( "mat*moveableTrans= %s")% mat.asString() );
    trans.translate(&fixedCenter);
    LOG(BF( "fixedTrans = %s")% trans.asString() );
    this->_Transform = trans*mat;
}


#if 0
/*!
  Return the rootMeanSquare difference between the two
  collections of points.
*/
double	SuperposeEngine_O::sumOfSquaresOfDifferences(ScorerState_sp scorerState )
{ _G();
    core::IntArray_O::iterator		itFixed;
    core::IntArray_O::iterator		ititMoved;
    Vector3				moved, diff;
    double				sum;

    ASSERTNOTNULL(this->_MoveableIndices);
    ASSERTNOTNULL(this->_FixedIndices);
    LOG_SCORE(scorerState,BF("SuperposeEngine_O::sumOfSquaresOfDifferences{"));
#if	DEBUG_SCORE_EVALUATION
    LOG_SCORE(scorerState,BF("There are %d moveable indices")%this->_MoveableIndices->size() );
    LOG_SCORE(scorerState,BF("There are %d fixed indices")%this->_FixedIndices->size() );
    for ( ititMoved = this->_MoveableIndices->begin(),
	      itFixed = this->_FixedIndices->begin();
	  itFixed != this->_FixedIndices->end();
	  ititMoved++, itFixed++ )
    {
	LOG_SCORE(BF(scorerState,"Superpose fixed%s with moveable%s")%
		  this->_FixedCoordinates->getElement(*itFixed).asString()
		  % this->_MoveableCoordinates->getElement(*ititMoved).asString()
	    );
    }
    LOG_SCORE(BF(scorerState,"SuperposeEngine_O calculated transform\n%s")% this->_Transform.asString());
#endif
    double result = this->sumOfSquaresOfDifferences();
    LOG_SCORE(scorerState,BF("SuperposeEngine_O::returning result=%lf")% result);
    LOG_SCORE(scorerState,BF("SuperposeEngine_O::sumOfSquaresOfDifferences}"));
    return result;
}
#endif

/*!
  Return the rootMeanSquare difference between the two
  collections of points.
*/
double	SuperposeEngine_O::sumOfSquaresOfDifferences()
{ _G();
    core::IntArray_O::iterator		itFixed;
    core::IntArray_O::iterator		ititMoved;
    Vector3				moved, diff;
    double				sum;

    LOG(BF("Dump of coordinates being that I will superpose") );
    ASSERTNOTNULL(this->_MoveableIndices);
    ASSERTNOTNULL(this->_FixedIndices);
#ifdef	DEBUG_ON
    LOG(BF("There are %d moveable indices") % this->_MoveableIndices->size()  );
    LOG(BF("There are %d fixed indices") % this->_FixedIndices->size()  );
    for ( ititMoved = this->_MoveableIndices->begin(),
	      itFixed = this->_FixedIndices->begin();
	  itFixed != this->_FixedIndices->end();
	  ititMoved++, itFixed++ )
    {
	LOG(BF("Superpose fixed%s with moveable%s")
	    % this->_FixedCoordinates->getElement(*itFixed).asString()
	    % this->_MoveableCoordinates->getElement(*ititMoved).asString() );
    }
#endif
    ASSERT(this->_FixedIndices->size()!=0);
    sum = 0.0;
    for ( ititMoved = this->_MoveableIndices->begin(),
	      itFixed = this->_FixedIndices->begin();
	  itFixed != this->_FixedIndices->end();
	  ititMoved++, itFixed++ ) {
	moved = this->_Transform.multiplyByVector3(this->_MoveableCoordinates->getElement(*ititMoved));
	diff = this->_FixedCoordinates->getElement(*itFixed)-moved;
	sum = sum + diff.dotProduct(diff);
    }
    LOG(BF("Calculated sum=%lf") % sum  );
    return sum;
}

//
//	rootMeanSquareDifference
//
//	Return the rootMeanSquare difference between the two
//	collections of points.
double	SuperposeEngine_O::rootMeanSquareDifference()
{ _G();
    double sumOfSquares = this->sumOfSquaresOfDifferences();
    LOG(BF("Number of moveable indices = %d") % this->_MoveableIndices->size() );
    ASSERT_gt(this->_MoveableIndices->size(),0);
    LOG(BF("sumOfSquares = %lf") % sumOfSquares );
    double meanSumOfSquares = sumOfSquares/this->_MoveableIndices->size();
    LOG(BF("meanSumOfSquares = %lf") % meanSumOfSquares );
    double rms = sqrt(meanSumOfSquares);
    LOG(BF("rms = %lf") % rms );
    return rms;
}






void	SuperposeEngine_O::setFixedPoints( core::IntArray_sp fi, geom::CoordinateArray_sp fc )
{ _G();
    LOG(BF("SuperposeEngine_O::setFixedPoints --> number of points=%d") % fc->size()  );
    if ( !((fi->size()>=3) && ((fi->size()<=fc->size()))) )
    {
	stringstream ss;
	ss << "In SuperposeEngine setFixedPoints - number of indices fi->size() = ";
	ss << fi->size() << " must be at least 3"<< std::endl;
	ss << "Number of coordinates --> fc->size() = ";
	ss << fc->size() << " must be greator than or equal to fi->size(";
	ss << fi->size() << ")" << std::endl;
	LOG(BF("%s") % ss.str().c_str()  );
	SIMPLE_ERROR(BF(ss.str()));
    }
    this->_FixedIndices = fi;
    this->_FixedCoordinates = fc;
}



void	SuperposeEngine_O::setFixedAllPoints( geom::CoordinateArray_sp fc )
{_G();
    core::IntArray_O::iterator	ia;
    uint		ii;
    this->_FixedCoordinates = fc;
    this->_FixedIndices = core::IntArray_O::create(fc->size(),_lisp);
    ii = 0;
    for ( ii=0; ii<fc->size(); ii++ )
    {
        this->_FixedIndices->put(ii,ii);
        LOG(BF("setFixedAllPoints index@%d ") % ii  );
    }
}



void	SuperposeEngine_O::setMoveablePoints( core::IntArray_sp mi, geom::CoordinateArray_sp mc )
{ _G();
    VectorVector3s::iterator	cur;
    LOG(BF("SuperposeEngine_O::setMoveablePoints --> number of points=%d") % mc->size()  );
    ASSERTP( (mi->size()>=3), "There must be at least three indices" );
    ASSERTP( mi->size()<=mc->size(), "There must be at least as many coordinates as indices");
    this->_MoveableIndices = mi;
    this->_MoveableCoordinates = mc->copy();
}

void	SuperposeEngine_O::setMoveableAllPoints( geom::CoordinateArray_sp mc )
{_G();
    core::IntArray_O::iterator	ia;
    uint		ii;
    ASSERTF(mc->size()>=3,BF("You must have at least three moveable points and there are only %d") % mc->size() );
    this->_MoveableCoordinates = mc;
    this->_MoveableIndices = core::IntArray_O::create(mc->size(),_lisp);
    for ( ii=0; ii<mc->size(); ii++ )
    {
        this->_MoveableIndices->put(ii,ii);
        LOG(BF("setMoveableAllPoints index@%d ") % ii  );
    }
}







Matrix	SuperposeEngine_O::superpose()
{_G();
    LOG(BF("SuperposeEngine_O::superpose()") );
    this->doSuperpose();
    LOG(BF("Carried out superpose and the transform is:%s") % (this->_Transform.asString().c_str() ) );
    return this->_Transform;
}

void SuperposeSelectedAtoms_O::exposeCando(core::Lisp_sp e)
{
    core::class_<SuperposeSelectedAtoms_O>()
	;

}

    void SuperposeSelectedAtoms_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,SuperposeSelectedAtoms,"","",_lisp)
	;
#endif
#endif //]
}






void SuperposeSelectedAtoms_O::initialize()
{_G();
    this->Base::initialize();
    this->_Matter = _Nil<Matter_O>();
    this->_SuperposeAtoms.clear();
}


#ifdef XML_ARCHIVE
void SuperposeSelectedAtoms_O::archiveBase(core::ArchiveP node)
{_G();
    this->Base::archiveBase(node);
    node->attribute("matter",this->_Matter);
    node->archiveVector0("superposeAtoms",this->_SuperposeAtoms);
}
#endif


void SuperposeSelectedAtoms_O::updateSuperposeAtoms()
{_G();
    Loop lAtoms;
    this->_SuperposeAtoms.clear();
    lAtoms.loopTopGoal(this->_Matter,ATOMS);
    uint numAtoms = 0;
    uint numSelectedAtoms = 0;
    SYMBOL_EXPORT_SC_(ChemPkg,superpose);
    core::Symbol_sp superposeSymbol = _sym_superpose;
    while ( lAtoms.advance() )
    {
	Atom_sp a = lAtoms.getAtom();
	if ( a->getPropertyOrDefault(superposeSymbol,_lisp->_true()).isTrue() )
	{
	    this->_SuperposeAtoms.push_back(a);
	    numSelectedAtoms++;
	}
	numAtoms++;
    }
    ASSERT_gt(numAtoms,0);
    ASSERT_gt(numSelectedAtoms,0);
    ASSERT_gt(this->_SuperposeAtoms.size(),0);
}


void SuperposeSelectedAtoms_O::setMatter(Matter_sp matter)
{_G();
    this->_Matter = matter->copy();
    this->_SuperposeAtoms.clear();
}

void SuperposeSelectedAtoms_O::setMatterWithSelectedAtoms(Matter_sp matter)
{_G();
    this->setMatter(matter);
    this->updateSuperposeAtoms();
}



geom::CoordinateArray_sp SuperposeSelectedAtoms_O::extractCoordinates(Matter_sp matter)
{_G();
    if ( !this->_Matter->equal(matter) )
    {
	SIMPLE_ERROR(BF("The Matters are not equal"));
    }
    this->_Matter->transferCoordinates(matter);
    ASSERT_gt(this->_SuperposeAtoms.size(),0);
    geom::CoordinateArray_sp coords = geom::CoordinateArray_O::create(this->_SuperposeAtoms.size());
    uint idx = 0;
    for (gctools::Vec0<Atom_sp>::iterator ai=this->_SuperposeAtoms.begin(); ai!=this->_SuperposeAtoms.end(); ai++, idx++ )
    {
	coords->setElement(idx,(*ai)->getPosition());
    }
    return coords;
}


Matter_sp SuperposeSelectedAtoms_O::getMatter()
{_G();
    return this->_Matter;
}

void SuperposeSelectedAtoms_O::copyMatterCoordinatesIntoFixedCoordinates(Matter_sp matter)
{_G();
    geom::CoordinateArray_sp ca = this->extractCoordinates(matter);
    ASSERT_gt(ca->size(),0);
    this->setFixedAllPoints(ca);
}

void SuperposeSelectedAtoms_O::copyMatterCoordinatesIntoMoveableCoordinates(Matter_sp matter)
{_G();
    geom::CoordinateArray_sp ca = this->extractCoordinates(matter);
    ASSERT_gt(ca->size(),0);
    this->setMoveableAllPoints(ca);
}













};

