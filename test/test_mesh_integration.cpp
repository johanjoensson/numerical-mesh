#include "numerical-mesh.h"
#include "numerical-mesh-integration.h"
#include <gtest/gtest.h>

#define TOL 1e-10

TEST(LinearTrapz, TestLinearTrapz1D)
{
    Linear_mesh<1, double> m(0., 5., 10);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = trapezoidal_integral<double>(m, 0., 10., 1., f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(LinearSimpson, TestLinearSimpson1D)
{
    Linear_mesh<1, double> m(0., 5., 10);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = simpson_integral<double>(m, 0., 10., 1., f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(LinearCorrected3, TestLinearCorrected31D)
{
    Linear_mesh<1, double> m(0., 5., 10);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = corrected_trapezoidal_integral<double, 3>(m, 0., 10., 1., f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(QuadraticTrapz, TestQuadraticTrapz1D)
{
    Quadratic_mesh<1, double> m(0., 5., 10);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = trapezoidal_integral<double>(m, 0., 10., 1., f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(QuadraticSimpson, TestQuadraticSimpson1D)
{
    Quadratic_mesh<1, double> m(0., 5., 10);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = simpson_integral<double>(m, 0., 10., 1., f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(QuadraticCorrected3, TestQuadraticCorrected31D)
{
    Quadratic_mesh<1, double> m(0., 5., 10);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = corrected_trapezoidal_integral<double, 3>(m, 0., 10., 1., f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(ExåponentialTrapz, TestExponentialTrapz1D)
{
    Exponential_mesh<1, double> m(0., 5., 0.000005, 10);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = trapezoidal_integral<double>(m, 0., 10., 1., f);
    ASSERT_NEAR(integral, 25., TOL) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(ExponentialSimpson, TestExponentialSimpson1D)
{
    Exponential_mesh<1, double> m(0., 5., 0.005, 200);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = simpson_integral<double>(m, 0., 200., 1., f);
    ASSERT_NEAR(integral, 25., TOL) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(ExponetialCorrected3, TestExponentialCorrected31D)
{
    Exponential_mesh<1, double> m(0., 5., 0.001, 10);
    auto f = [] (const double& x) { return 5 + x*0;};
    auto integral = corrected_trapezoidal_integral<double, 3>(m, 0., 10., 1., f);
    ASSERT_NEAR(integral, 25., TOL) << "Integral 0 -> 5 of 5 dr != 25";
}
