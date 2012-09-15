
/******************************************************************************

  blas.hpp, v 1.0.2 2011.10.14

  Basic linear algebra subprograms. 

  Copyright (c) 2011
  Dzmitry Hlindzich <dzmitry.hlindzich@ziti.uni-heidelberg.de>
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1.  Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
  2.  Redistributions in binary form must reproduce the above copyright
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

// Suppress boost::numeric::ublas C4127 warning under MSVC.
#ifdef _MSC_VER
#   pragma warning(push)
#   pragma warning(disable:4127)
#endif // _MSC_VER
#   include <boost/numeric/ublas/matrix.hpp>
#   include <boost/numeric/ublas/lu.hpp>
#ifdef _MSC_VER
#   pragma warning(pop)
#endif // _MSC_VER

#include "bo/extended_math.hpp"

using namespace boost::numeric::ublas;

namespace bo {
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
    size_t res = lu_factorize(A, pm);
    if (res != 0)
        return false;

    // Create identity matrix of "inverse"
    inverse.assign(identity_matrix<T> (A.size1()));

    // Backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
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
            pm_sign *= -1;
        }
    }

    return pm_sign;
}

// Calculate the determinant of the input matrix.
template<class T>
double determinant(const matrix<T>& input)
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


// Eigenvector decomposition for real symmetric matrices.
// Returns a vector of the eigenvalues, sorted in non-decreasing order.
// The corresponding eigenvectors are stored in the columns of the
// matrix A.
template <class T>
std::vector<T> eigen_analysis(matrix<T>& A)
{
    unsigned int n = A.size1();

    std::vector<T> d(n);

    // Initialize the vector.
    for (int j = 0; j < n; ++j)
    {
        d[j] = A(n - 1, j);
    }

    BOOST_ASSERT(A.size2() == n);

    // Check for a square matrix.
    if (A.size2() != n)
    {
        return d;
    }

    std::vector<T> e(n, T(0));

    // Iterating.
    for (int i = n - 1; i > 0; --i) 
    {
        T scale(0);

        for (int k = 0; k < i; ++k)
        {
            scale += std::fabs(d[k]);
        }

        if (scale == T(0))
        {            
            e[i] = d[i - 1];

            for (int j = 0; j < i; ++j)
            {
                d[j] = A(i - 1, j); 
                A(i, j) = A(j, i) = T(0); 
            }

            d[i] = T(0);
        } 
        else
        {
            T h(0);
            T invscale = T(1.0 / scale);

            for (int k = 0; k < i; ++k)
            {
                d[k] *= invscale;
                h += square(d[k]);
            }

            T f = d[i - 1];
            T g = (f > 0) ? -std::sqrt(h) : std::sqrt(h);
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;

            for (int j = 0; j < i; ++j)
            {
                e[j] = T(0);
            }

            for (int j = 0; j < i; ++j) 
            {
                f = d[j];
                A(j, i) = f;
                g = e[j] + f * A(j, j);

                for (int k = j+1; k < i; k++)
                {
                    g += A(k, j) * d[k];
                    e[k] += A(k, j) * f;
                }

                e[j] = g;
            }

            f = T(0);
            T invh = T(1.0 / h);

            for (int j = 0; j < i; ++j)
            {
                e[j] *= invh;
                f += e[j] * d[j];
            }

            T hh = f / (h + h);

            for (int j = 0; j < i; j++)
            {
                e[j] -= hh * d[j];
            }

            for (int j = 0; j < i; j++)
            {
                f = d[j];
                g = e[j];

                for (int k = j; k < i; k++)
                {
                    A(k, j) -= f * e[k] + g * d[k];
                }

                d[j] = A(i - 1, j);
                A(i, j) = T(0);
            }

            d[i] = h;
        }
    }

    // Doing some magic.
    for (int i = 0; i < n - 1; ++i)
    {
        A(n - 1, i) = A(i, i);
        A(i, i) = 1;

        T h = d[i+1];

        if (h != T(0))
        {
            T invh = T(1.0 / h);

            for (int k = 0; k <= i; ++k)
            {
                d[k] = A(k, i + 1) * invh;
            }

            for (int j = 0; j <= i; ++j)
            {
                T g(0);

                for (int k = 0; k <= i; ++k)
                {
                    g += A(k, i + 1) * A(k, j);
                }

                for (int k = 0; k <= i; ++k)
                {
                    A(k, j) -= g * d[k];
                }

            }
        }

        for (int k = 0; k <= i; ++k)
        {
            A(k, i + 1) = T(0);
        }    
    }

    for (int j = 0; j < n; ++j) 
    {
        d[j] = A(n - 1, j);
        A(n - 1, j) = T(0);
    }

    A(n - 1, n - 1) = 1;

    // QL.
    for (int i = 1; i < n; ++i)
    {
        e[i - 1] = e[i];
    }

    e[n - 1] = T(0);

    T f(0), tmp(0);

    const T eps = T(std::pow(T(2), -52));

    for (int l = 0; l < n; l++) 
    {
        tmp = std::max(tmp, std::fabs(d[l]) + std::fabs(e[l]));
        int m = l;

        while (m < n)
        {
            if (std::fabs(e[m]) <= eps * tmp)
                break;
            ++m;
        }

        if (m > l)
        {
            do
            {
                T g = d[l];
                T p = (d[l + 1] - g) / (e[l] + e[l]);
                T r = T(std::sqrt(square(p) + T(1)));

                if (p < T(0))
                {
                    r = -r;
                }

                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                T dl1 = d[l + 1];
                T h = g - d[l];

                for (int i = l + 2; i < n; ++i)
                {
                    d[i] -= h;
                }

                f += h;
                p = d[m];   
                T c(1), c2(1), c3(1);

                T el1 = e[l + 1];
                T s(0), s2(0);

                for (int i = m - 1; i >= l; --i)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;

                    r = T(std::sqrt(square(p) + square(e[i])));

                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);

                    for (int k = 0; k < n; ++k)
                    {
                        h = A(k, i + 1);
                        A(k, i + 1) = s * A(k, i) + c * h;
                        A(k, i) = c * A(k, i) - s * h;
                    }
                }

                p = -s * s2 * c3 * el1 * e[l] / dl1;

                e[l] = s * p;
                d[l] = c * p;

            }
            while (std::fabs(e[l]) > eps * tmp);
        }

        d[l] += f;
        e[l] = T(0);
    }

    // Sort
    for (int i = 0; i < n - 1; ++i)
    {
        int k = i;
        T p = d[i];

        for (int j = i + 1; j < n; ++j)
        {
            if (d[j] < p)
            {
                k = j;
                p = d[j];
            }
        }

        if (k == i) continue;

        d[k] = d[i];
        d[i] = p;

        for (int j = 0; j < n; ++j)
        {
            p = A(j, i);;
            A(j, i) = A(j, k);
            A(j, k) = p;
        }
    }
}

} // namespace blas
} // namespace bo

#endif // BLAS_HPP_F74A6974_6444_4C40_BFE7_75ADEC15B7E6_
