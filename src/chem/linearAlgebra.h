       
       

//
// (C) 2004 Christian E. Schafmeister
//



#ifndef LinearAlgebra_H
#define	LinearAlgebra_H

#include <string>
#include "energyFunction.h"
#include "largeSquareMatrix.h"

namespace       chem {



	void unconventionalModifiedCholeskySymbolicFactorization(chem::AbstractLargeSquareMatrix_sp m,
							chem::AbstractLargeSquareMatrix_sp l);
	void unconventionalModifiedCholeskyFactorization(chem::SparseLargeSquareMatrix_sp m,
				chem::SparseLargeSquareMatrix_sp l, chem::NVector_sp tempVector );
	void backSubstituteLDLt(chem::AbstractLargeSquareMatrix_sp ldlt,
				chem::NVector_sp x, chem::NVector_sp b );


};

#endif
