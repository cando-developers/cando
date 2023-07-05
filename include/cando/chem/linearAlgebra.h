/*
    File: linearAlgebra.h
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



#ifndef LinearAlgebra_H
#define	LinearAlgebra_H

#include <string>
#include <cando/chem/energyFunction.h>
#include <cando/chem/largeSquareMatrix.h>

namespace       chem_old {
	void unconventionalModifiedCholeskySymbolicFactorization(chem::AbstractLargeSquareMatrix_sp m,
							chem::AbstractLargeSquareMatrix_sp l);
	void unconventionalModifiedCholeskyFactorization(chem::SparseLargeSquareMatrix_sp m,
				chem::SparseLargeSquareMatrix_sp l, chem::NVector_sp tempVector );
	void backSubstituteLDLt(chem::AbstractLargeSquareMatrix_sp ldlt,
				chem::NVector_sp x, chem::NVector_sp b );
};

namespace       chem {
	void unconventionalModifiedCholeskySymbolicFactorization(chem::AbstractLargeSquareMatrix_sp m,
							chem::AbstractLargeSquareMatrix_sp l);
	void unconventionalModifiedCholeskyFactorization(chem::SparseLargeSquareMatrix_sp m,
				chem::SparseLargeSquareMatrix_sp l, chem::NVector_sp tempVector );
	void backSubstituteLDLt(chem::AbstractLargeSquareMatrix_sp ldlt,
				chem::NVector_sp x, chem::NVector_sp b );
};

#endif
