#include <numerical-mesh-integration.h>
#include <gtest/gtest.h>

#define TOL 5e-15

TEST(LinearTrapz, TestLinearTrapz1D)
{
    Linear_mesh<1, double> m(1., 6., 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = trapezoidal_integral<double>(m, f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(LinearSimpson, TestLinearSimpson1D)
{
    Linear_mesh<1, double> m(1., 6., 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = simpson_integral<double>(m, f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(LinearCorrected3, TestLinearCorrected31D)
{
    Linear_mesh<1, double> m(1., 6., 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = corrected_trapezoidal_integral<double, 3>(m, f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(QuadraticTrapz, TestQuadraticTrapz1D)
{
    Quadratic_mesh<1, double> m(1., 6., 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = trapezoidal_integral<double>(m, f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(QuadraticSimpson, TestQuadraticSimpson1D)
{
    Quadratic_mesh<1, double> m(1., 6., 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = simpson_integral<double>(m, f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(QuadraticCorrected3, TestQuadraticCorrected31D)
{
    Quadratic_mesh<1, double> m(1., 6., 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = corrected_trapezoidal_integral<double, 3>(m, f);
    ASSERT_DOUBLE_EQ(integral, 25.) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(ExponentialTrapz, TestExponentialTrapz1D)
{
    Exponential_mesh<1, double> m(1., 6., 1e-8, 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = trapezoidal_integral<double>(m, f);
    ASSERT_NEAR(integral, 25., TOL) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(ExponentialSimpson, TestExponentialSimpson1D)
{
    Exponential_mesh<1, double> m(1., 6., 1e-8, 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = simpson_integral<double>(m, f);
    ASSERT_NEAR(integral, 25., TOL) << "Integral 0 -> 5 of 5 dr != 25";
}

TEST(ExponetialCorrected3, TestExponentialCorrected31D)
{
    Exponential_mesh<1, double> m(1., 6., 1e-8, 11);
    auto f = [] (const Mesh_base<1, double>::mesh_point& p) { return 5 + p.r()*0;};
    auto integral = corrected_trapezoidal_integral<double, 3>(m, f);
    ASSERT_NEAR(integral, 25., TOL) << "Integral 0 -> 5 of 5 dr != 25";
}
