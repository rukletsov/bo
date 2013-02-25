
/******************************************************************************

  distances_3d.hpp, v 1.0.6 2013.02.25

  Methods and algorithms for calculating distances between varios objects in
  a 3D-space.

  Copyright (c) 2011 - 2013
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

#ifndef DISTANCES_3D_HPP_B5895686_0C10_449A_9DB3_03BEDBE065FB_
#define DISTANCES_3D_HPP_B5895686_0C10_449A_9DB3_03BEDBE065FB_

#include "bo/vector.hpp"
#include "bo/extended_math.hpp"

namespace bo {
namespace methods {

namespace detail {

// Functions for processing different regions of the st-plane for the global minimum
// of the point-to-face distance function Q(s, t). These functions are used only in
// Mesh::distance_to_face() function. For more information see comments on this function.

// The simplest case: the global minimum of Q(s, t) is located inside the triangle.
// Just perform a postponed division.
inline
void region0_helper(double& s, double& t, const double& denominator)
{
    double inv_factor = 1. / denominator;
    s *= inv_factor;
    t *= inv_factor;
}

// The global minimum of Q(s, t) is on the other side of the (E0-E1)-edge of the
// triangle. Therefore, the closest point lies on the (E0-E1)-edge and can be found
// by minimizing Q(s, 1-s) function. It has the only minimum at (c+e-b-d)/(a-2b+c).
// Since (a-2b+c) is always greater than 0 for a triangle, consider only (c+e-b-d)
// and according to its value determine the s-parameter of the closest point.
template <typename T>
void region1_helper(double& s, double& t, const T& a, const T& b, const T& c,
                    const T& d, const T& e)
{
    T numerator = c + e - b - d;

    if (numerator <= T(0))
        s = 0.;
    else
    {
        T denominator = a - T(2) * b + c;
        s = ((numerator >= denominator) ? 1. :
                                          (double(numerator) / double(denominator)));
    }

    t = 1. - s;
}

// The global minimum of Q(s, t) is on the other side of the (E1-B)-edge of the
// triangle. Therefore, the closest point lies on the (E1-B)-edge and can be found by
// minimizing Q(0, t) function. It has the only minimum at -e/c. Since c is always
// greater than 0 for a triangle, consider only -e and according to its value determine
// the t-parameter of the closest point.
template <typename T>
void region3_helper(double& s, double& t, const T& c, const T& e)
{
    s = 0.;

    if (e >= T(0))
        t = 0.;
    else
        t = ((-e >= c) ? 1. : (double(-e) / double(c)));
}

// The global minimum of Q(s, t) is on the other side of the (E0-B)-edge of the
// triangle. Therefore, the closest point lies on the (E0-B)-edge and can be found by
// minimizing Q(s, 0) function. It has the only minimum at -d/a. Since a is always
// greater than 0 for a triangle, consider only -d and according to its value determine
// the s-parameter of the closest point.
template <typename T>
void region5_helper(double& s, double& t, const T& a, const T& d)
{
    t = 0.;

    if (d >= T(0))
        s = 0.;
    else
        s = ((-d >= a) ? 1. : (double(-d) / double(a)));
}

// The global minimum of Q(s, t) is between the prolongations of (E1-B)-edge and
// (E0-E1)-edge. However the closest point can lie on one of the both edges and in
// point (0, 1). The negative gradient of Q(s, t) cannot point inside the triangle.
// Consider two vectors which represent edges: (1, -1) for (E0-E1)-edge and (0, -1)
// for (E1-B)-edge. The dot products of the positive gradient Grad(Q(s, t)) with
// these vectors shows the direction of the -Grad(Q(0, 1)) which implies what edge is
// closer to the minimum. If [(1, -1) dot Grad(Q(0, 1) < 0] then the closest point
// lies on the (E0-E1)-edge, and formula from region1 can be used to determine it;
// if [(0, -1) dot Grad(Q(0, 1)) < 0] then the closest point lies on the (E1-B)-edge
// and formula from region3 can be used; otherwise the closest point is (0, 1). Two
// dot products cannot be less than 0 simultaneously.
template <typename T>
void region2_helper(double& s, double& t, const T& a, const T& b, const T& c,
                    const T& d, const T& e)
{
    T tmp0 = b + d;
    T tmp1 = c + e;
    if (tmp1 > tmp0)
    {
        // Minimum on the (E0-E1)-edge. Same as region1.
        T numerator = tmp1 - tmp0;
        T denominator = a - T(2) * b + c;
        s = ((numerator >= denominator) ? 1. :
                                          (double(numerator) / double(denominator)));
        t = 1. - s;
    }
    else
    {
        // Minimum on the (E1-B)-edge. Same as region3.
        s = 0.;
        if (tmp1 <= T(0))
            // Minimum in (0, 1)
            t = 1.;
        else
            // Minimum somewhere on the (E1-B)-edge.
            t = ((e >= T(0)) ? 0. : (double(-e) / double(c)));
    }
}

// The global minimum of Q(s, t) is between the prolongations of (E0-B)-edge and
// (E0-E1)-edge. However the closest point can lie on one of the both edges and in
// point (1, 0). The negative gradient of Q(s, t) cannot point inside the triangle.
// Consider two vectors which represent edges: (-1, 1) for (E0-E1)-edge and (-1, 0)
// for (E0-B)-edge. The dot products of the positive gradient Grad(Q(s, t)) with
// these vectors shows the direction of the -Grad(Q(1, 0)) which implies what edge is
// closer to the minimum. If [(-1, 1) dot Grad(Q(1, 0) < 0] then the closest point
// lies on the (E0-E1)-edge, and formula from region1 can be used to determine it;
// if [(-1, 0) dot Grad(Q(1, 0)) < 0] then the closest point lies on the (E0-B)-edge
// and formula from region5 can be used; otherwise the closest point is (1, 0). Two
// dot products cannot be less than 0 simultaneously.
template <typename T>
void region6_helper(double& s, double& t, const T& a, const T& b, const T& c,
                    const T& d, const T& e)
{
    T tmp0 = b + e;
    T tmp1 = a + d;
    if (tmp1 > tmp0)
    {
        // Minimum on the (E0-E1)-edge. Similar to region1, but Q(1-t, t) is used.
        T numerator = tmp1 - tmp0;
        T denominator = a - T(2) * b + c;
        t = ((numerator >= denominator) ? 1. :
                                          (double(numerator) / double(denominator)));
        s = 1. - s;
    }
    else
    {
        // Minimum on the (E0-B)-edge. Same as region5.
        t = 0.;
        if (tmp1 <= T(0))
            // Minimum in (1, 0)
            s = 1.;
        else
            // Minimum somewhere on the (E0-B)-edge.
            s = ((d >= T(0)) ? 0. : (double(-d) / double(a)));
    }
}

// The global minimum of Q(s, t) is between the prolongations of (E1-B)-edge and
// (E0-B)-edge. However the closest point can lie on one of the both edges and in
// point (0, 0). The negative gradient of Q(s, t) cannot point inside the triangle.
// Consider two vectors which represent edges: (1, 0) for (E0-B)-edge and (0, 1) for
// (E1-B)-edge. The dot products of the positive gradient Grad(Q(s, t)) with these
// vectors shows the direction of the -Grad(Q(0, 0)) which implies what edge is
// closer to the minimum. If [(1, 0) dot Grad(Q(0, 0) < 0] then the closest point
// lies on the (E0-B)-edge, and formula from region5 can be used to determine it;
// if [(0, 1) dot Grad(Q(0, 0)) < 0] then the closest point lies on the (E1-B)-edge
// and formula from region3 can be used; otherwise the closest point is (0, 0). Two
// dot products cannot be less than 0 simultaneously.
template <typename T>
void region4_helper(double& s, double& t, const T& a, const T& c, const T& d, const T& e)
{
    if (d < T(0))
        // Minimum on the (E0-B)-edge. Same as region5.
        s = ((-d >= a) ? 1. : (double(-d) / double(a)));
    else
        // Minimum on the (E1-B)-edge including (0, 0).
        s = 0.;

    if (e < T(0))
        // Minimum on the (E1-B)-edge. Same as region3.
        t = ((-e >= c) ? 1. : (double(-e) / double(c)));
    else
        // Minimum on the (E0-B)-edge including (0, 0).
        t = 0.;
}

} // namespace detail


// Computes the euclidean distance between two given points in arbitrary space.
// Points may be of integer type, however the result is always double.
template <typename T, std::size_t N>
double euclidean_distance_d(const bo::Vector<T, N>& point1,
                            const bo::Vector<T, N>& point2)
{
    double distance = (point1 - point2).eucl_norm();
    return distance;
}

// Computes the euclidean distance between two given points in arbitrary space.
// Points are supposed to be of some real type, which is supported by std::sqrt
// function. Otherwise ambiguous call or undefined function for RealType may happen.
template <typename RealType, std::size_t N>
RealType euclidean_distance(const bo::Vector<RealType, N>& point1,
                            const bo::Vector<RealType, N>& point2)
{
    bo::Vector<RealType, N> point_diff = (point1 - point2);
    RealType distance = std::sqrt(point_diff * point_diff);
    return distance;
}

// Computes the taxicab (or L1) distance between two given points in arbitrary space.
// Points may be of any type T, for which bo::Vector<T, N>::taxicab_norm() works.
template <typename T, std::size_t N>
T taxicab_distance(const bo::Vector<T, N>& point1, const bo::Vector<T, N>& point2)
{
    T distance = (point1 - point2).taxicab_norm();
    return distance;
}

// Computes Chebyshev (or L inf) distance between two given points in arbitrary space.
// Points may be of any type T, for which bo::Vector<T, N>::maximum_norm() works.
template <typename T, std::size_t N>
T chebyshev_distance(const bo::Vector<T, N>& point1, const bo::Vector<T, N>& point2)
{
    T distance = (point1 - point2).maximum_norm();
    return distance;
}


// Finds the closest point on the given triangle to the given point. The type of the
// vertex is supposed to be one of the floating point types (float, double, long double).
// Integer types can be also used but due to rounding the result won't be precise.
// However, a custom type T fully compatible to double can be used without any errors.
// The requirements are: common arithmetic operations support, constructor T(double)
// exists and (preferably) does no rounding. For more information see
//    http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
// or check available documents in the "/papers" directory.
template <typename T>
bo::Vector<T, 3> find_closest_point_on_triangle(const bo::Vector<T, 3>& P,
    const bo::Vector<T, 3>& corner1, const bo::Vector<T, 3>& corner2,
    const bo::Vector<T, 3>& corner3)
{
    // Prepare parametrized representation of the triangle
    // T(s, t) = B + sE0 + tE1.
    const bo::Vector<T, 3> B = corner1;
    const bo::Vector<T, 3> E0 = corner2 - B;
    const bo::Vector<T, 3> E1 = corner3 - B;

    // Prepare coefficients for distance function Q(s, t).
    T a = E0 * E0;
    T b = E0 * E1;
    T c = E1 * E1;
    T d = E0 * (B - P);
    T e = E1 * (B - P);
    // T f = (B - P) * (B - P);

    // Compute the global mimimum of the distance function Q(s, t).
    double denominator = a * c - b * b;
    double s = b * e - c * d;
    double t = b * d - a * e;

    // Determine one of the 7 regions, where the global minimum is located.
    if (s + t <= denominator)
    {
        // Global minimum is on the "left" side of the st-plane, including triangle.
        if (s < 0.)
        {
            if (t < 0.)
                detail::region4_helper(s, t, a, c, d, e);
            else
                detail::region3_helper(s, t, c, e);
        }
        else if (t < 0.)
            detail::region5_helper(s, t, a, d);
        else
            detail::region0_helper(s, t, denominator);
    }
    else
    {
        // Global minimum is on the "right" side of the st-plane, not including triangle.
        if (s < 0.)
            detail::region2_helper(s, t, a, b, c, d, e);
        else if (t < 0.)
            detail::region6_helper(s, t, a, b, c, d, e);
        else
            detail::region1_helper(s, t, a, b, c, d, e);
    }

    // After processing the corresponding region, in s and t we have the closest point
    // on the triangle to the given point.
    bo::Vector<T, 3> closest_point = B + T(s) * E0 + T(t) * E1;
    return closest_point;
}


// Computes the projection of a point q onto a plane given by a point and a normal.
template <typename T>
bo::Vector<T, 3> project_point_onto_plane(const bo::Vector<T, 3>& point,
    const bo::Vector<T, 3>& plane_origin, const bo::Vector<T, 3>& plane_norm)
{
    // Ensure the plane normal is normalized.
    bo::Vector<T, 3> unit_norm(plane_norm.normalized(), 3);

    // Find the distance to plane along the normal.
    T dist = unit_norm * (point - plane_origin);

    // Multiply the unit normal by the distance, and subtract this vector from the point.
    bo::Vector<T, 3> retvalue = point - (unit_norm * dist);

    return retvalue;
}

} // namespace methods
} // namespace bo

#endif // DISTANCES_3D_HPP_B5895686_0C10_449A_9DB3_03BEDBE065FB_
