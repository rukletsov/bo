
/******************************************************************************

    blas.hpp, v 1.0.0 2011.03.16

    Basic linear algebra subprograms. 

    Copyright (c) 2011
    Dzmitry Hlindzich <dzmitry.hlindzich@ziti.uni-heidelberg.de>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    1.	Redistributions of source code must retain the above copyright
	    notice, this list of conditions and the following disclaimer.
    2.	Redistributions in binary form must reproduce the above copyright
	    notice, this list of conditions and the following disclaimer in the
	    documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
    OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
    SUCH DAMAGE.

*******************************************************************************/

#ifndef BLAS_HPP_F74A6974_6444_4C40_BFE7_75ADEC15B7E6_
#define BLAS_HPP_F74A6974_6444_4C40_BFE7_75ADEC15B7E6_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric::ublas;

namespace blas {

// Matrix inversion routine.
// Use lu_factorize and lu_substitute in uBLAS to invert a matrix. 
template<class T>
bool invert_matrix(const matrix<T>& input, matrix<T>& inverse)
{
	// Create a working copy of the input
	matrix<T> A(input);

	// Create a permutation matrix for the LU-factorization
	permutation_matrix<std::size_t> pm(A.size1());

	// Perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// Create identity matrix of "inverse"
	inverse.assign(identity_matrix<T> (A.size1()));

	// Backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

// Calculate the determinant of the input matrix.
template<class T>
double determinant(const matrix<T>& input )
{
	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	permutation_matrix<std::size_t> pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);

	double det=1.0;

	if (res != 0 )
	{
		det = 0.0;
	} 
	else
	{
		// multiply by elements on diagonal
		for(int i = 0; i < A.size1(); ++i) 
			det *= A(i,i); 
		det = det * determinant_sign( pm );
	}

	return det;
}

// Define the sign of the determinant using the given permutation matrix.
inline
int determinant_sign(const permutation_matrix<std::size_t>& pm)
{
	int pm_sign=1;
	std::size_t size = pm.size();
	for (std::size_t i = 0; i < size; ++i)
	{
		if (i != pm(i)) 
		{
			// swap_rows would swap a pair of rows here, so we change sign
			pm_sign *= -1.0;
		}
	}

	return pm_sign;
}

} // namespace blas

#endif // BLAS_HPP_F74A6974_6444_4C40_BFE7_75ADEC15B7E6_
