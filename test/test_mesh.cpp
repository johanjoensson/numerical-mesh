#include "numerical-mesh.h"
#include <gtest/gtest.h>

#define TOL 1e-10

TEST(LinearMesh, TestLinear1D)
{
    Linear_mesh<1, double> m(1., 6., 11);
    ASSERT_DOUBLE_EQ(m.r(10), 6.) << "r(10) = 0.5*10 + 1 != 6.";
    ASSERT_DOUBLE_EQ(m.r(-10), -4.) << "r(-10) = 0.5*(-10) + 1 != -4.";
    ASSERT_DOUBLE_EQ(m.r2(10), 36.) << "r^2(10) = (0.5*10 + 1)^2 != 36.";
    ASSERT_DOUBLE_EQ(m.dr(5), 0.5) << "d/dx(0.5*x + 1)|x=5 != 0.5";
    ASSERT_DOUBLE_EQ(m.dr(-5), 0.5) << "d/dx(0.5*x + 1)|x=-5 != 0.5";
    ASSERT_DOUBLE_EQ(m(10), 6.) << "m(10) = 0.5*10 + 1 != 6.";
    ASSERT_DOUBLE_EQ(m(-10), -4.) << "m(-10) = 0.5*(-10) + 1 != -4.";

}

TEST(QuadraticMesh, TestQuadratic1D)
{
    Quadratic_mesh<1, double> m(1., 6., 11);
    ASSERT_DOUBLE_EQ(m.r(10), 6.) << "r(10) = 0.05*10^2 + 1 != 6.";
    ASSERT_DOUBLE_EQ(m.r(-10), -4.) << "r(-10) = -0.05*10^2 + 1 != -4.";
    ASSERT_DOUBLE_EQ(m.r2(10), 36.) << "r^2(10) = (0.05*10^2 + 1)^2 != 36.";
    ASSERT_DOUBLE_EQ(m.dr(5), 0.5) << "d/dx(0.05*x^2 + 1)|x=5 != 0.5";
    ASSERT_DOUBLE_EQ(m.dr(-5), 0.5) << "-d/dx(0.05*x^2 + 1)|x=-5 != 0.5";
    ASSERT_DOUBLE_EQ(m(10), 6.) << "m(10) = 0.05*10^2 + 1 != 6.";
    ASSERT_DOUBLE_EQ(m(-10), -4.) << "m(-10) = -0.05*10^2 + 1 != -4.";

}

TEST(ExponentialMesh, TestExponential1D)
{
    Exponential_mesh<1, double> m(1., 6., 0.02, 11);
    ASSERT_DOUBLE_EQ(m.r(10), 6.) << "r(10) = 0.02*(e^(0.02*10) - 1) + 1 != 6.";
    ASSERT_DOUBLE_EQ(m.r(-10), -4.) << "r(-10) = -0.02*(e^(0.02*10) - 1) + 1 != -4.";
    ASSERT_NEAR(m.r2(10), 36., TOL) << "r^2(10) = (0.02*(e^(0.02*10) - 1) + 1)^2 != 36.";
    ASSERT_NEAR(m.dr(5), 0.49916763786, TOL) << "d/dx(22.58..*(e^(0.02x) - 1) + 1)|x=5 != 0.49916763786";
    ASSERT_NEAR(m.dr(-5), 0.49916763786, TOL) << "-d/dx(22.58..*(e^(0.02x) - 1) + 1)|x=-5 != 0.49916763786";
    ASSERT_DOUBLE_EQ(m(10), 6.) << "m(10) = 0.02*(e^(0.02*10) - 1) + 1 != 6.";
    ASSERT_DOUBLE_EQ(m(-10), -4.) << "m(-10) = -0.02*(e^(0.02*10) - 1) + 1 != -4.";
}

TEST(LinearMesh, TestLinear5D)
{
    Linear_mesh<5, double> m({1., 1., 1., 1., 1.}, {6., 6., 6., 6., 6.}, {11, 11, 11, 11, 11});

    std::ranges::for_each(m.r({10., 10., 10., 10., 10.}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, 6.) << "r(10) = 0.5*10 + 1 != 6.";
        });
    std::ranges::for_each(m.r({-10, -10, -10, -10, -10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, -4) << "r(-10) = 0.5*(-10) + 1 != -4.";
        });
    std::ranges::for_each(m.r2({10, 10, 10, 10, 10}),
        [] (const double r2i)
        {
            ASSERT_DOUBLE_EQ(r2i, 36.) << "r^2(10) = (0.5*10 + 1)^2 != 36.";
        });
    std::ranges::for_each(m.dr({5, 5, 5, 5, 5}),
        [] (const double dri)
        {
            ASSERT_DOUBLE_EQ(dri, 0.5) << "d/dx(0.5*x + 1)|x=5 != 0.5";
        });
    std::ranges::for_each(m.dr({-5, -5, -5, -5, -5}),
        [] (const double dri)
        {
            ASSERT_DOUBLE_EQ(dri, 0.5) << "d/dx(0.5*x + 1)|x=-5 != 0.5";
        });
    std::ranges::for_each(m({10, 10, 10, 10, 10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, 6.) << "r(10) = 0.5*10 + 1 != 6.";
        });
    std::ranges::for_each(m({-10, -10, -10, -10, -10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, -4.) << "r(-10) = 0.5*(-10) + 1 != -4.";
        });

}

TEST(QuadraticMesh, TestQuadratic5D)
{
    Quadratic_mesh<5, double> m({1., 1., 1., 1., 1.}, {6., 6., 6., 6., 6.}, {11, 11, 11, 11, 11});

    std::ranges::for_each(m.r({10., 10., 10., 10., 10.}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, 6.) << "r(10) = 0.5*10 + 1 != 6.";
        });
    std::ranges::for_each(m.r({-10, -10, -10, -10, -10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, -4) << "r(-10) = 0.5*(-10) + 1 != -4.";
        });
    std::ranges::for_each(m.r2({10, 10, 10, 10, 10}),
        [] (const double r2i)
        {
            ASSERT_DOUBLE_EQ(r2i, 36.) << "r^2(10) = (0.5*10 + 1)^2 != 36.";
        });
    std::ranges::for_each(m.dr({5, 5, 5, 5, 5}),
        [] (const double dri)
        {
            ASSERT_DOUBLE_EQ(dri, 0.5) << "d/dx(0.5*x + 1)|x=5 != 0.5";
        });
    std::ranges::for_each(m.dr({-5, -5, -5, -5, -5}),
        [] (const double dri)
        {
            ASSERT_DOUBLE_EQ(dri, 0.5) << "d/dx(0.5*x + 1)|x=-5 != 0.5";
        });
    std::ranges::for_each(m({10, 10, 10, 10, 10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, 6.) << "r(10) = 0.5*10 + 1 != 6.";
        });
    std::ranges::for_each(m({-10, -10, -10, -10, -10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, -4.) << "r(-10) = 0.5*(-10) + 1 != -4.";
        });
}

TEST(ExponentialMesh, TestExponential5D)
{
    Exponential_mesh<5, double> m({1., 1., 1., 1., 1.}, {6., 6., 6., 6., 6.}, {0.02, 0.02, 0.02, 0.02, 0.02}, {11, 11, 11, 11, 11});

    std::ranges::for_each(m.r({10., 10., 10., 10., 10.}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, 6.) << "r(10) = 0.5*10 + 1 != 6.";
        });
    std::ranges::for_each(m.r({-10, -10, -10, -10, -10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, -4) << "r(-10) = 0.5*(-10) + 1 != -4.";
        });
    std::ranges::for_each(m.r2({10, 10, 10, 10, 10}),
        [] (const double r2i)
        {
            ASSERT_DOUBLE_EQ(r2i, 36.) << "r^2(10) = (0.5*10 + 1)^2 != 36.";
        });
    std::ranges::for_each(m.dr({5, 5, 5, 5, 5}),
        [] (const double dri)
        {
            ASSERT_NEAR(dri, 0.49916763786, TOL) << "d/dx(22.58..*(e^(0.02x) - 1) + 1)|x=5 != 0.49916763786";
        });
    std::ranges::for_each(m.dr({-5, -5, -5, -5, -5}),
        [] (const double dri)
        {
            ASSERT_NEAR(dri, 0.49916763786, TOL) << "-d/dx(22.58..*(e^(0.02x) - 1) + 1)|x=-5 != 0.49916763786";
        });
    std::ranges::for_each(m({10, 10, 10, 10, 10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, 6.) << "r(10) = 0.5*10 + 1 != 6.";
        });
    std::ranges::for_each(m({-10, -10, -10, -10, -10}),
        [] (const double ri)
        {
            ASSERT_DOUBLE_EQ(ri, -4.) << "r(-10) = 0.5*(-10) + 1 != -4.";
        });
}

TEST(MeshTest, TestEvaluate1D)
{
    auto func = [] (const double r)
    {
        return std::sin(2*M_PI*r);
    };
    Linear_mesh<1, double> ml(0, 1, 101);
    Quadratic_mesh<1, double> mq(0, 1, 101);
    Exponential_mesh<1, double> me(0, 1, 0.02, 101);
    auto linear_res = ml.evaluate(func);
    auto quadratic_res = mq.evaluate(func);
    auto exponential_res = me.evaluate(func);

    ASSERT_NEAR(linear_res.front(), 0, TOL);
    ASSERT_NEAR(quadratic_res.front(), 0, TOL);
    ASSERT_NEAR(exponential_res.front(), 0, TOL);

    ASSERT_NEAR(linear_res.back(), 0, TOL);
    ASSERT_NEAR(quadratic_res.back(), 0, TOL);
    ASSERT_NEAR(exponential_res.back(), 0, TOL);

}
