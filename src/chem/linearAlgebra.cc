/*
    File: linearAlgebra.cc
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
       
#define	TURN_DEBUG_OFF
#include <clasp/core/common.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/linearAlgebra.h>


namespace chem {



void	unconventionalModifiedCholeskySymbolicFactorization(
			AbstractLargeSquareMatrix_sp 	m,
			AbstractLargeSquareMatrix_sp	ldlt )
{
int	dim, j;
#define	TAO	10.0	// See Schlick paper page 10, top paragraph

    dim = m->dimension();
		//
		// Calculate each column and the diagonal
		//
    for ( j=0; j<dim; j++ ) { 
	ldlt->insertElement(j,j);	// Diagonal element is always there
	for ( int i=j+1; i<dim; i++ ) {	// i is X-coordinate
	    if ( m->hasElement(i,j) ) {
		ldlt->insertElement(i,j);
	    };
	    if ( j>0 && ldlt->hasElement(i,j-1) ) {
		ldlt->insertElement(i,j);
	    }
	}
    }
    ldlt->insertionIsComplete();
}

void	unconventionalModifiedCholeskyFactorization(
			SparseLargeSquareMatrix_sp 	m,
			SparseLargeSquareMatrix_sp	ldlt,
			NVector_sp			nvKSum)
{
uint	dim;
uint	k;
uint	j;
double	dk, ksum, ljk, thetaSquared;
double	epsilon, betaSquared, delta,cjk;
uint	cjkIndex;
double	djbar, djtilde, dj, theta, cik, cij;
double	ddim;
int	rowFirstIndex, rowLastIndex, ii,xII;
double	ldltJk;
#define	TAO	10.0	// See Schlick paper page 10, top paragraph

//_lisp->profiler().timer(core::timerPreconditioner).start();
//_lisp->profiler().timer(core::timerPreconditionerSolver).start();
//_lisp->profiler().timer(core::timerPreconditionerSolverFactor).start();
    dim = m->dimension();
    ddim = (double)dim;
		//
		// Calculate each column and the diagonal
		//
    epsilon = m->maxAbsValue();
    LOG("epsilon = %lf" , epsilon  );
    betaSquared = epsilon/(sqrt(ddim*(ddim-1.0)));
    delta = 1e-9;
    LOG("betaSquared = %lf" , betaSquared  );
    for ( j=0; j<dim; j++ ) { 
	ksum = 0;
	for ( k=0; k<j; k++ ) {	// k is a Y coordinate
	    cjkIndex = ldlt->indexFromCoordinatesOrUndefinedUnsignedInt(j,k);
	    if ( cjkIndex == UndefinedUnsignedInt ) continue;
	    cjk = ldlt->getAtIndex(cjkIndex);
	    dk = ldlt->element(k,k);
	    ASSERTP(dk!=0.0, "unconventionalModifiedCholeskyFactorization>> "
				"about to divide by zero");
	    ljk = cjk/dk;
	    ldlt->setElement(j,k,ljk);
	    LOG("Storing l[%d',%d']=%lf" , j+1 , k+1 , ljk );
	    ksum += ljk*cjk;
	}
	djbar = m->element(j,j)-ksum;
		//
		// Calculate theta for the coeff's and theta for the diagonal
		//
	theta = 0.0;
//
// The slow stuff is from here.
//
//        _lisp->profiler().timer(core::timerPreconditionerSolverFactorFocus).start();
	LOG("Starting ksum calculation" );
#if 1 //[    Faster implementation, for now comment out
	nvKSum->zero();	// initialize ksum vector to zero
		// k is Y coordinate
	for (k=0; k<j; k++ ) 
	{ 
	    rowFirstIndex = ldlt->indexOfFirstElementAtOrAfterX(j,k);
	    LOG("Check if there is a ldlt(j,k) element" );
			// If this row doesn't have an element at x=j
			// then ldlt[j,k] = 0 and there is nothing to
			// accumulate
	    if ( ldlt->columnForIndex(rowFirstIndex) > j ) continue;
	    LOG("Getting rowLastIndex" );
	    rowLastIndex = ldlt->indexOfLastElementOnRow(k);
	    LOG("Getting ldltJk" );
	    ldltJk = ldlt->getAtIndex(rowFirstIndex);
	    LOG("Accumulating nvKSum" );
				//ii indexes along X
	    for (ii=rowFirstIndex; ii<= rowLastIndex; ii++ ) {
		cik = ldlt->getAtIndex(ii);
		xII = ldlt->columnForIndex(ii);
		nvKSum->addToElement(xII,cik*ldltJk);
	    }
	}
	if ( j+1<dim ) 
	{
	    LOG("Setting up to calculate cij j+1(%d) columns(%d)" , j+1 , ldlt->columns() );
	    rowFirstIndex = ldlt->indexOfFirstElementAtOrAfterX(j+1,j);
	    rowLastIndex = ldlt->indexOfLastElementOnRow(j);
			    // ii indexes along X-coordinate
	    LOG("Calculating cij" );
	    for ( ii=rowFirstIndex; ii<=rowLastIndex; ii++ ) 
	    { 
		xII = ldlt->columnForIndex(ii);
		LOG("xII=%d" , xII );
		cij = m->element(xII,j)-nvKSum->element(xII);
		LOG("cij=%lf" , cij );
		ldlt->setElement(xII,j,cij);
		LOG("put cij into ldlt at element %d,%d" , xII , j );
		theta = MAX(theta,fabs(cij));
	    }
	}
#else //][
	for ( i=j+1; i<dim; i++ ) 
	{	// i is X-coordinate
	    ksum = 0.0;
	    for ( k=0; k<j; k++ ) 
	    {	// k is Y coordinate
		cik = ldlt->element(i,k); // Temporarily store cik in l
		if ( cik == 0.0 ) continue;	// MEISTER OPTIMIZE HERE
		ksum += ldlt->element(j,k)*cik;
	    }
	    cij = m->element(i,j)-ksum;
	    hasij = ldlt->hasElement(i,j);
	    if ( !hasij && cij!= 0.0 ) {
		stringstream ss;
		ss << "Want to write "<<cij<<" to ldlt["<<i<<","<<j<<"] but it isn't in ldlt";
		SIMPLE_ERROR(("%s") , ss.str().c_str());
	    }
	    if ( !hasij ) continue;
	    LOG("Storing c[%d',%d']=%lf" , i+1 , j+1 , cij  );
	    ldlt->setElement(i,j,cij);
	    theta = MAX(theta,fabs(cij));
	}
#endif //]
//        _lisp->profiler().timer(core::timerPreconditionerSolverFactorFocus).stop();
//
// To here
//
	LOG("Calculated theta=%lf" , theta  );
	thetaSquared = theta*theta;
	double tsdbs = thetaSquared/betaSquared;
	djtilde = djbar + TAO;
	if ( djtilde > delta ) {
	    if ( tsdbs == 0 ) {
		dj = djtilde;
	    } else {
		dj = MIN(djtilde,thetaSquared/betaSquared);
	    }
	} else if ( djtilde < -delta ) {
	    if ( tsdbs == 0 ) {
		dj = djtilde;
	    } else {
	        dj = MAX(djtilde,-thetaSquared/betaSquared);
	    }
	} else {
	    dj = delta;
	}
	LOG("Storing d[%d']=%lf" , j+1 , dj );
	ldlt->setElement(j,j,dj);
    }
//    _lisp->profiler().timer(core::timerPreconditionerSolverFactor).stop();
//    _lisp->profiler().timer(core::timerPreconditionerSolver).stop();
//    _lisp->profiler().timer(core::timerPreconditioner).stop();
}

SYMBOL_EXPORT_SC_(ChemPkg,trap_back_substitute_ldl);

void	backSubstituteLDLt(AbstractLargeSquareMatrix_sp ldlt,
                           NVector_sp s1, NVector_sp b )
{
  ASSERT(ldlt->_Triangle==SymmetricDiagonalLower);
  NVector_sp	tx;
  int	        x,y;
  double	sum, dd;
  tx = NVector_O::create(b->size());
  if (gc::IsA<SparseLargeSquareMatrix_sp>(ldlt)) {
    SparseLargeSquareMatrix_sp slsm = gc::As_unsafe<SparseLargeSquareMatrix_sp>(ldlt);
      // Lower diagonal
    for ( y=0; y<b->size(); y++ ) {
      sum = 0.0;
      slsm->loopOverRowLowerDiagonal(y,0,y,[&s1,&sum](uint col, uint row, double value) {
        double element = s1->element(col);
        sum += value * element;
      } );
      s1->setElement(y,b->element(y)-sum);
    }
    for ( x=0; x<s1->size(); x++ ) {
      dd = ldlt->element(x,x);
      ASSERTP(dd!=0.0,"backSubstituteLDLt diagonal element shouldn't be zero");
      tx->setElement(x,s1->element(x)/dd);
    }
    for ( y=b->size()-1; y>=0; y-- ) {
      sum = 0.0;
      slsm->loopOverColumnLowerDiagonal(y,y+1,b->size(),[&s1,&sum](uint col, uint row, double value) {
        double element = s1->element(col);
        sum += value * element;
      } );
      s1->setElement(y,tx->element(y)-sum);
    }
  } else {
    for ( y=0; y<b->size(); y++ ) {
      sum = 0.0;
      for ( x=0; x<y; x++ ) {
        sum += ldlt->element(x,y)*s1->element(x);
      }
      s1->setElement(y,b->element(y)-sum);
    }
    for ( x=0; x<s1->size(); x++ ) {
      dd = ldlt->element(x,x);
      ASSERTP(dd!=0.0,"backSubstituteLDLt diagonal element shouldn't be zero");
      tx->setElement(x,s1->element(x)/dd);
    }
    for ( y=b->size()-1; y >= 0; y-- ) {
      sum = 0.0;
      for ( x=y+1; x<b->size(); x++ ) {
        double val = ldlt->element(x,y);
        double element = s1->element(x);
        sum += val*element;
      }
      s1->setElement(y,tx->element(y)-sum);
    }
  }
}


};
