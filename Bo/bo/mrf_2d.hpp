
/******************************************************************************

  mrf.hpp, v 0.0.1 2012.09.06

  Markov random field model for regular 2D lattice.

  Copyright (c) 2009 - 2012
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

#ifndef MRF_HPP_1140ED81_1E3E_4AEB_AFBB_6CA70FE3EF9B_
#define MRF_HPP_1140ED81_1E3E_4AEB_AFBB_6CA70FE3EF9B_

#include <cmath>
#include <functional>
#include <boost/random.hpp>

#include "bo/raw_image_2d.hpp"

namespace bo {

template <typename NodeType, typename DataType, typename RealType>
class MRF2D
{
public:
    static const int kError = -1;

public:
    typedef RawImage2D<NodeType> RandomLattice;
    typedef RawImage2D<DataType> DataLattice;
    typedef std::binary_function<DataType, NodeType, RealType> LikelihoodEnergy;
    typedef std::binary_function<NodeType, NodeType, RealType> PriorEnergy;

//    typedef boost::variate_generator<boost::mt19937&, boost::uniform_real<double> >
//        Generator;

    MRF2D(const RandomLattice& initial_configuration, const DataLattice& observation,
          LikelihoodEnergy likelihood, PriorEnergy prior);

    const Image& current() const;

    double next_icm_iteration();
    double next_mmd_iteration(Generator& rng, const double& temperature,
                              bool is_modified, const double& mmd_probab);

    double p() const;
    double compute_full_energy() const;
    double compute_local_energy(ValType new_val, size_t row, size_t col);

    static ValType diff(const Image& image, co::Index ind1, co::Index ind2);
    static double pow2(ValType val);
    static ValType abs(ValType val);
    static double fi_gm(ValType val);
    static double fi_gr(ValType val);
    static double tau_threshold(ValType val, double tau);

protected:
    double compute_alpha_(const Image& model) const;
    double compute_beta_(const Image& model) const;

    double right_clique_(const Image& model, size_t row, size_t col) const;
    double down_clique_(const Image& model, size_t row, size_t col) const;
    double left_clique_(const Image& model, size_t row, size_t col) const;
    double up_clique_(const Image& model, size_t row, size_t col) const;

    double diff_clique_(const Image& model, const Image& initial,
                        size_t row, size_t col) const;

    double clique_fun_(const Image& image, co::Index ind1, co::Index ind2) const;

    double neighb_penalty(const Image& model, co::Index ind1, co::Index ind2) const;
    double edge_penalty(const Image& model, co::Index ind1, co::Index ind2) const;

    Pixels get_values() const;
    Pixels get_neighbour_values(const Image& image, size_t row, size_t col) const;

private:
    Image observation_;
    Image current_;

    // Disallow copy and construct
    MRF2D(const MRF2D&);
    MRF2D& operator=(const MRF2D&);
};


template <typename ValType>
GibbsModel<ValType>::GibbsModel(const Image& observation, const Image& model,
                                const double& alpha, const double& beta,
                                const double& gamma, int colours,
                                const double& thres_border):
    observation_(observation),
    current_(model),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma),
    colours_(colours),
    thres_border_(thres_border)
{ }

template <typename ValType> inline
typename const GibbsModel<ValType>::Image& GibbsModel<ValType>::current() const
{
    return current_;
}

// Implements an iteration of the ICM algorithm
template <typename ValType>
double GibbsModel<ValType>::next_icm_iteration()
{
    // Clone a previous state.
    Image new_state = current_;

    // Determine a set of intensities to use in the ICM algorithm.
    Pixels values = get_values();

    // Do an inner iteration pixelwise.
    for (size_t row = 0; row < current_.size(); ++row)
    {
        for (size_t col = 0; col < current_[row].size(); ++col)
        {   // Search for a value which gives the minimum of energy.
            double min_energy = std::numeric_limits<double>::max();
            int min_value_index = -1;
            //Pixels values = get_neighbour_values(current_, row, col);

            for (size_t i = 0; i < values.size(); ++i)
            {
                double energy = compute_local_energy(values[i], row, col);
                if (energy < min_energy)
                {
                    min_energy = energy;
                    min_value_index = static_cast<int>(i);
                }
            }

            // Apply the best value.
            new_state[row][col] = values[min_value_index];
        }
    }

    // Compute an improved difference and proceed to a new state.
    double old_p = p();
    current_ = new_state;
    double new_p = p();

    return
        (new_p - old_p);
}

// Implements an iteration of the MMD algorithm.
template <typename ValType>
double GibbsModel<ValType>::next_mmd_iteration(Generator& rng, const double& temperature,
                                               bool is_modified, const double& mmd_probab)
{
    // Clone a previous state.
    Image new_state = current_;

    // If a modified version of MD algo (MMD) is used, then the decision probability is
    // fixed and obtained as a parameter. Otherwise, it is generated every time.
    double decision_probab = log(mmd_probab);

    // Current temperature. See algorithm description for details.
    double t = temperature;

    // Determine a set of intensities (possible values) for each pixel.
    Pixels values = get_values();

    // Do an inner iteration pixelwise.
    for (size_t row = 0; row < current_.size(); ++row)
    {
        for (size_t col = 0; col < current_[row].size(); ++col)
        {   // Randomly generate a new state for the current pixel.
            size_t index = static_cast<size_t>(co::round(rng() * (values.size() - 1)));
            ValType new_val = values[index];

            // if classical version is used, generate the decision probability.
            if (!is_modified)
                decision_probab = log(rng());

            // Compute local energy and accept it or reject.
            if (decision_probab <=
                    (compute_local_energy(current_[row][col], row, col) -
                    compute_local_energy(new_val, row, col)) / t
               )
            {   // Accept new state.
                new_state[row][col] = new_val;
            }
        }
    }

    // Compute an improved difference and proceed to a new state.
    double old_energy = compute_full_energy();
    current_ = new_state;
    double new_energy = compute_full_energy();

    return
        (new_energy - old_energy);
}

// Returns a probability of a current MRF state.
template <typename ValType> inline
double GibbsModel<ValType>::p() const
{
    return
        exp(static_cast<double>(0.0 - compute_full_energy()));
}

// Computes energy of a current MRF state.
template <typename ValType>
double GibbsModel<ValType>::compute_full_energy() const
{
    double alphagamma_item = compute_alpha_(current_);
    double beta_item = compute_beta_(current_);
    double energy = alphagamma_item + beta_item;

    return energy;
}

// Computes an energy of a state with one pixel [row, col] intensity changed to "new_val".
// Border checks for pixel's border belonging are done inside necessary clique functions.
template <typename ValType>
double GibbsModel<ValType>::compute_local_energy(ValType new_val, size_t row, size_t col)
{
    double energy = 0.0;

    // IMPORTANT: in order to minimize computational time by avoiding copying
    // we temporary modify the current MRF state (current_). Before the function
    // returns the state is restored. This can lead to memory corruption in case
    // of concurrent usage. To make this function thread-safe, make it const
    // (either copy data or rewrite clique functions so they support values).
    ValType old_val = current_[row][col];
    current_[row][col] = new_val;

    // Compute a difference of pixel intensity of model and observation.
    energy += diff_clique_(current_, observation_, row, col);

    // Compute a difference of neighbour pixel intensities. Border overrun is checeked
    // inside clique functions.
    energy += right_clique_(current_, row, col);
    energy += down_clique_(current_, row, col);
    energy += left_clique_(current_, row, col);
    energy += up_clique_(current_, row, col);

    // Restore the state.
    current_[row][col] = old_val;

    return energy;
}

// This function computes first item of energy which goes with coefficient "alpha".
// It takes the sum over all first order cliques (every two adjacent pixels),
// performing subtraction of pixel intensities and then squaring the result.
// In order to iterate over all first order cliques we can for every pixel
// consider RIGHT and DOWN neighbour except for border right and border down ones.
// For down-right pixel we should do nothing. Since clique functions check border
// overrun we can simply compute right and down cliques for every pixel. For the
// undefined results zero energy is returned.
template <typename ValType>
double GibbsModel<ValType>::compute_alpha_(const Image& model) const
{
    double sum = 0.0;

    // Process rows excluding the last one. Exclusion is done inside the clique functions.
    for (size_t row = 0; row < model.size(); ++row)
    {
        // Process a row till (last-1) pixel. See the note above.
        for (size_t col = 0; col < model[row].size(); ++col)
        {
            sum += right_clique_(model, row, col);
            sum += down_clique_(model, row, col);
        }
    }

    return sum;
}

// This function computes second item of energy which goes with coefficient "beta".
// It squares the difference between pixel intensity on input image and on current
// MRF model.
template <typename ValType>
double GibbsModel<ValType>::compute_beta_(const Image& model) const
{
    // Check if the dimentions of observed image and MRF model are the same.
    if ( (model.size() != observation_.size()) ||
         (model[0].size() != observation_[0].size()) )
    {
        co::errprint(kAppName, "Input image and MRF model are of different sizes.");
        return kError;
    }

    double sum = 0.0;

    // Process pixel by pixel
    for (size_t row = 0; row < observation_.size(); ++row)
        for (size_t col = 0; col < observation_[row].size(); ++col)
            sum += diff_clique_(model, observation_, row, col);

    return sum;
}

// Next five functions compute neighbouring cliques. First four (right, down, left, up)
// count the difference of intensities of given and a neighbour pixel of an image.
// Fifth one counts the difference of intensities of the same pixel on two different
// images. It is supposed, that one image is the current state of MRF model and the other
// one is an initial observation. All these functions do check border overrun.
template <typename ValType> inline
double GibbsModel<ValType>::right_clique_(const Image& model, size_t row, size_t col) const
{
    if (col != model[row].size() - 1)
    {
        return
            clique_fun_(model, std::make_pair(row, col), std::make_pair(row, col + 1));
    }
    else
        return 0.0;
}

template <typename ValType> inline
double GibbsModel<ValType>::down_clique_(const Image& model, size_t row, size_t col) const
{
    if (row != model.size() - 1)
    {
        return
            clique_fun_(model, std::make_pair(row, col), std::make_pair(row + 1, col));
    }
    else
        return 0.0;
}

template <typename ValType> inline
double GibbsModel<ValType>::left_clique_(const Image& model, size_t row, size_t col) const
{
    if (col != 0)
    {
        return
            clique_fun_(model, std::make_pair(row, col), std::make_pair(row, col - 1));
    }
    else
        return 0.0;
}

template <typename ValType> inline
double GibbsModel<ValType>::up_clique_(const Image& model, size_t row, size_t col) const
{
    if (row != 0)
    {
        return
            clique_fun_(model, std::make_pair(row, col), std::make_pair(row - 1, col));
    }
    else
        return 0.0;
}

template <typename ValType> inline
double GibbsModel<ValType>::diff_clique_(const Image& model, const Image& initial,
                                         size_t row, size_t col) const
{
    // Either
    //    pow2(val) or
    //    abs(val) or
    //    fi_gm(val).
    //return
    //    fi_gm(initial[row][col] - model[row][col]);

    double retvalue = abs(initial[row][col] - model[row][col]);

    return
        (beta_ * retvalue);
}

template <typename ValType> inline
double GibbsModel<ValType>::clique_fun_(const Image& image, co::Index ind1, co::Index ind2) const
{
    return
        (alpha_ * neighb_penalty(image, ind1, ind2) +
         gamma_ * edge_penalty(image, ind1, ind2));
}

// Next function represents mathematical functionals applied to couplings in
// the GibbsModel.
template <typename ValType> inline
double GibbsModel<ValType>::neighb_penalty(const Image& model, co::Index ind1, co::Index ind2) const
{
    return
        pow2(diff(model, ind1, ind2));
}

template <typename ValType> inline
double GibbsModel<ValType>::edge_penalty(const Image& model, co::Index ind1, co::Index ind2) const
{
    return
        tau_threshold(diff(model, ind1, ind2), thres_border_);
}

// Determines a set of possible intensities for the underlying image.
template <typename ValType>
typename GibbsModel<ValType>::Pixels GibbsModel<ValType>::get_values() const
{
    Pixels retvalue;

    int values_quantity = colours_;
    int max_val = values_quantity - 1;

    retvalue.reserve(values_quantity);
    for (int i = 0; i < values_quantity; ++i)
        retvalue.push_back(static_cast<ValType>(i) / static_cast<ValType>(max_val));

    return retvalue;
}

template <typename ValType> inline
ValType GibbsModel<ValType>::diff(const Image& image, co::Index ind1, co::Index ind2)
{
    return
        (image[ind1.first][ind1.second] - image[ind2.first][ind2.second]);
}

template <typename ValType> inline
double GibbsModel<ValType>::pow2(ValType val)
{
    return
        pow(static_cast<double>(val), 2);
}

template <typename ValType> inline
ValType GibbsModel<ValType>::abs(ValType val)
{
    return
        ::abs(static_cast<double>(val));
}

// Returns "(x^2) / (1+x^2)".
template <typename ValType> inline
double GibbsModel<ValType>::fi_gm(ValType val)
{
    return
        (static_cast<double>(val * val) / static_cast<double>(1.0 + val * val));
}

// Returns "2ln(cosh(x))".
template <typename ValType> inline
double GibbsModel<ValType>::fi_gr(ValType val)
{
    return
        2.0 * log(
            (exp(static_cast<double>(val)) + exp(static_cast<double>(0.0 - val))) / 2.0);
}

// Returns "|x|*(1-psi(s,t)) + tau * psi(s,t)",
// where "psi(s,t) = 1, if |x| > tau,
//                   0, if |x| <= tau".
template <typename ValType> inline
double GibbsModel<ValType>::tau_threshold(ValType val, double tau)
{
    return std::min(abs(val), tau);
}

} // namespace bo

#endif // MRF_HPP_1140ED81_1E3E_4AEB_AFBB_6CA70FE3EF9B_
