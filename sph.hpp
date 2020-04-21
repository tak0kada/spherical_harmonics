#pragma once

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace sph
{

    using real_t = double;

    // see the document for the definition
    inline real_t sph_harm(const unsigned int& l, const int& m, const real_t& theta, const real_t& phi) {
        using namespace boost::math;
        constexpr real_t sqrt2{double_constants::root_two};

        if (m > 0) {
            return sqrt2 * (m % 2 == 0? 1: -1) * spherical_harmonic_r(l, m, theta, phi);
        }
        else if (m < 0) {
            return sqrt2 * (m % 2 == 0? 1: -1) * spherical_harmonic_i(l, -m, theta, phi);
        }
        else {// m == 0 {
            return spherical_harmonic_r(l, 0, theta, phi);
        }
    }

}
