#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/detail/tolerance_manip.hpp>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <array>
#include <cmath>
#include "sph.hpp"

inline double x(const double& theta, const double& phi)
{
    return std::sin(theta) * std::cos(phi);
}

inline double y(const double& theta, const double& phi)
{
    return std::sin(theta) * std::sin(phi);
}

inline double z(const double& theta, const double& phi)
{
    return std::cos(theta);
}


BOOST_AUTO_TEST_CASE(test_sph)
{
    using sph::sph_harm;
    constexpr double pi = boost::math::constants::pi<double>();

    constexpr std::size_t n_theta{50};
    constexpr std::size_t n_phi{100};

    std::array<double, n_theta> theta{};
    for (std::size_t i = 0; i < theta.size(); ++i)
    {
        theta[i] = pi * i / static_cast<double>(theta.size());
    }
    std::array<double, n_phi> phi{};
    for (std::size_t i = 0; i < phi.size(); ++i)
    {
        phi[i] = 2 * pi * i / static_cast<double>(phi.size());
    }

    #define CHECK_SPH(L, M, ans, tolerance) { \
        for (std::size_t i = 0; i < n_theta; ++i) { \
            for (std::size_t j = 0; j < n_phi; ++j) { \
                BOOST_TEST(sph::sph_harm(L, M, theta[i], phi[j]) == ans(theta[i], phi[j]), tolerance); \
            } \
        } \
    }

    {
        // l = 0, m = 0
        auto ans = [](double theta, double phi) -> double
        {
            return 0.5 / std::sqrt(pi);
        };
        CHECK_SPH(0, 0, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 1, m = -1
        auto ans = [](double& theta, double& phi) -> double
        {
            return std::sqrt(0.75 / pi) * y(theta, phi);
        };
        CHECK_SPH(1, -1, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 1, m = 0
        auto ans = [](double theta, double phi) -> double
        {
            return std::sqrt(0.75 / pi) * z(theta, phi);
        };
        CHECK_SPH(1, 0, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 1, m = 1
        auto ans = [](double theta, double phi) -> double
        {
            return std::sqrt(0.75 / pi) * x(theta, phi);
        };
        CHECK_SPH(1, 1, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 2, m = -2
        auto ans = [](double theta, double phi) -> double
        {
            return 0.5 * std::sqrt(15 / pi) * x(theta, phi) * y(theta, phi);
        };
        CHECK_SPH(2, -2, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 2, m = -1
        auto ans = [](double theta, double phi) -> double
        {
            return 0.5 * std::sqrt(15 / pi) * y(theta, phi) * z(theta, phi);
        };
        CHECK_SPH(2, -1, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 2, m = 0
        auto ans = [](double theta, double phi) -> double
        {
            const double x2 = std::pow(x(theta, phi), 2);
            const double y2 = std::pow(y(theta, phi), 2);
            const double z2 = std::pow(z(theta, phi), 2);
            return 0.25 * std::sqrt(5 / pi) * (-x2 - y2 + 2 * z2);
        };
        CHECK_SPH(2, 0, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 2, m = 1
        auto ans = [](double theta, double phi) -> double
        {
            return 0.5 * std::sqrt(15 / pi) * z(theta, phi) * x(theta, phi);
        };
        CHECK_SPH(2, 1, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 2, m = 2
        auto ans = [](double theta, double phi) -> double
        {
            const double x2 = std::pow(x(theta, phi), 2);
            const double y2 = std::pow(y(theta, phi), 2);
            return 0.25 * std::sqrt(15 / pi) * (x2 - y2);
        };
        CHECK_SPH(2, 2, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 3, m = -3
        auto ans = [](double theta, double phi) -> double
        {
            const double x2 = std::pow(x(theta, phi), 2);
            const double y2 = std::pow(y(theta, phi), 2);
            return 0.25 * std::sqrt(17.5 / pi) * (3 * x2 - y2) * y(theta, phi);
        };
        CHECK_SPH(3, -3, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 3, m = -2
        auto ans = [](double theta, double phi) -> double
        {
            return 0.5 * std::sqrt(105 / pi) * x(theta, phi) * y(theta, phi) * z(theta, phi);
        };
        CHECK_SPH(3, -2, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 3, m = -1
        auto ans = [](double theta, double phi) -> double
        {
            const double x2 = std::pow(x(theta, phi), 2);
            const double y2 = std::pow(y(theta, phi), 2);
            const double z2 = std::pow(z(theta, phi), 2);
            return 0.25 * std::sqrt(10.5 / pi) * y(theta, phi) * (4 * z2 - x2 - y2);
        };
        CHECK_SPH(3, -1, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 3, m = 0
        auto ans = [](double theta, double phi) -> double
        {
            const double x2 = std::pow(x(theta, phi), 2);
            const double y2 = std::pow(y(theta, phi), 2);
            const double z2 = std::pow(z(theta, phi), 2);
            return 0.25 * std::sqrt(7 / pi) * z(theta, phi) * (2 * z2 - 3 * x2 - 3 * y2);
        };
        CHECK_SPH(3, 0, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 3, m = 1
        auto ans = [](double theta, double phi) -> double
        {
            const double x2 = std::pow(x(theta, phi), 2);
            const double y2 = std::pow(y(theta, phi), 2);
            const double z2 = std::pow(z(theta, phi), 2);
            return 0.25 * std::sqrt(10.5 / pi) * x(theta, phi) * (4 * z2 - x2 - y2);
        };
        CHECK_SPH(3, 1, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 3, m = 2
        auto ans = [](double theta, double phi) -> double
        {
            const double x2 = std::pow(x(theta, phi), 2);
            const double y2 = std::pow(y(theta, phi), 2);
            return 0.25 * std::sqrt(105 / pi) * (x2 - y2) * z(theta, phi);
        };
        CHECK_SPH(3, 2, ans, boost::test_tools::tolerance(1e-12));
    }
    {
        // l = 3, m = 3
        auto ans = [](double theta, double phi) -> double
        {
            const double x2 = std::pow(x(theta, phi), 2);
            const double y2 = std::pow(y(theta, phi), 2);
            return 0.25 * std::sqrt(17.5 / pi) * (x2 - 3 * y2) * x(theta, phi);
        };
        CHECK_SPH(3, 3, ans, boost::test_tools::tolerance(1e-12));
    }
}
