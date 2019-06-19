## Associated Legendre polynomial
Associated Ledgendre polynomials are defined in terms of derivatives of Legendre polynomial,
$$
P\_{l}^{m}(x) := (-1)^{m} (1-x^{2})^{m/2} \frac{d^{m}}{dx^{m}} P\_{l}(x).
$$

Legendre polynomial again defined in terms of derivatives is,
$$
P\_{l}(x) := \frac{1}{2^{l} l!} \frac{d^{l}}{dx^{l}}\left[(x^{2}-1)^{l}\right].
$$

* Reference
* [wikipedia/Associated_Legendre_polynomials](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials)

## Spherical harmonics
Using the associated Legendre polynomials defined above, general definition of spherical harmonics is given as
$$
Y\_{l}^{m}(\theta, \phi ) := \sqrt{{(2l+1)\over 4\pi}{(l-m)!\over (l+m)!}} P\_{l}^{m}(\cos{\theta}) e^{i m \phi}.
$$

Then, Real form basis of spherical harmonics we are to use is
$$
Y\_{l}^{m} :=
  \begin{cases}
    & \sqrt{2} (-1)^{m} \text{Im}[Y\_{l}^{|m|}] & \text{if}\ m<0 \\\\
    & Y\_{l}^{0} & \text{if}\ m=0 \\\\
    & \sqrt{2} (-1)^{m} \text{Re}[Y\_{l}^{m}] & \text{if}\ m>0.
  \end{cases}
$$

* Reference
* [wikipedia/Spherical_harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics)

## Implementation
```cpp
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace sph
{

    using real_t = double;

    real_t sph_harm(const unsigned int& l, const int& m, const real_t& theta, const real_t& phi) {
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
```