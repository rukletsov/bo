
/******************************************************************************

  pca.hpp, v 1.0.0 2012.12.12

  Principal component analysis implementation using the covariance method.

  Copyright (c) 2012
  Alexander Rukletsov <rukletsov@gmail.com>
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

#ifndef PCA_HPP_C2EA373C_F360_43BA_BACB_4B26B78759BE_
#define PCA_HPP_C2EA373C_F360_43BA_BACB_4B26B78759BE_

#include <vector>
#include <algorithm>
#include <functional>
#include <boost/tuple/tuple.hpp>
#include <boost/array.hpp>
#include <boost/assert.hpp>

#include "bo/vector.hpp"
#include "bo/extended_math.hpp"
#include "bo/blas/blas.hpp"
#include "bo/blas/conversions.hpp"

namespace bo {
namespace blas {

using namespace boost::numeric::ublas;

// Performs PCA for the given data. Returns a tuple of eigenvalues and eigenvectors.
// Implemented as a functor in order to store induced typedefs inside.
template <typename RealType, std::size_t Dim>
struct PCA
{
    typedef Vector<RealType, Dim> Sample;
    typedef std::vector<Sample> Samples;

    typedef RealType EigenValue;
    typedef boost::array<EigenValue, Dim> EigenValues;
    typedef Vector<RealType, Dim> EigenVector;
    typedef boost::array<EigenVector, Dim> EigenVectors;
    typedef boost::tuples::tuple<EigenValues, EigenVectors> Result;

    Result operator() (Samples data)
    {
        typedef blas::bounded_vector<RealType, Dim> BlasVector;
        typedef std::vector<RealType> StdVector;
        typedef blas::bounded_matrix<RealType, Dim, Dim> Matrix;

        // Find the mean among the neighbours.
        Sample mean_value = bo::mean(data);

        // Calculate the deviations from mean.
        std::transform(data.begin(), data.end(), data.begin(),
                       std::bind2nd(std::minus<Sample>(), mean_value));

        // Initialize covariance matrix.
        Matrix covar_matrix = blas::zero_matrix<RealType>(Dim);

        // Iteratively compute the Dim x Dim covariance matrix (sample by sample).
        for (Samples::const_iterator pt = data.begin(); pt != data.end(); ++pt)
        {
            BlasVector vect = blas::from_bo_vector(*pt);
            covar_matrix += blas::outer_prod(vect, vect);
        }

        covar_matrix /= data.size();

        // Get the eigenvectors of the covariance matrix.
        StdVector covar_eigenvalues = blas::eigen_symmetric(covar_matrix);
        BOOST_ASSERT((covar_eigenvalues.size() == Dim) &&
                     "Eigenvalues count differs from samples' dimensions.");

        // Convert eigenvectors to the output format.
        EigenValues eigenvalues;
        std::copy(covar_eigenvalues.begin(), covar_eigenvalues.end(), eigenvalues.begin());

        // Extract eignvectors from the columns of the updated covariance matrix and
        // convert them to the output format.
        EigenVectors eigenvectors;
        for (std::size_t col_idx = 0; col_idx < Dim; ++col_idx)
            eigenvectors[col_idx] = blas::to_bo_vector(BlasVector(blas::column(covar_matrix, col_idx)));

        // Put eigenvalues and eigenvectors into a tuple and return.
        return
            boost::tuples::make_tuple(eigenvalues, eigenvectors);
    }
};

} // namespace blas
} // namespace bo

#endif PCA_HPP_C2EA373C_F360_43BA_BACB_4B26B78759BE_
