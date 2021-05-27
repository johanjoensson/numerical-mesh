#include <numerical-mesh.h>
#include <gtest/gtest.h>

#define TOL 5e-15

TEST(LinearMesh, TestLinear1D)
{
    Linear_mesh<1, double> m(1., 6., 11);
    ASSERT_DOUBLE_EQ(m.r(10), 6.) << "r(10) = 0.5*10 + 1 != 6.";
    ASSERT_DOUBLE_EQ(m.r(-10), -4.) << "r(-10) = 0.5*(-10) + 1 != -4.";
    ASSERT_DOUBLE_EQ(m.r2(10), 36.) << "r^2(10) = (0.5*10 + 1)^2 != 36.";
    ASSERT_DOUBLE_EQ(m.dr(5), 0.5) << "d/dx(0.5*x + 1)|x=5 != 0.5";
    ASSERT_DOUBLE_EQ(m.dr(-5), 0.5) << "d/dx(0.5*x + 1)|x=-5 != 0.5";
}

TEST(QuadraticMesh, TestQuadratic1D)
{
    Quadratic_mesh<1, double> m(1., 6., 11);
    ASSERT_DOUBLE_EQ(m.r(10), 6.) << "r(10) = 0.05*10^2 + 1 != 6.";
    ASSERT_DOUBLE_EQ(m.r(-10), -4.) << "r(-10) = -0.05*10^2 + 1 != -4.";
    ASSERT_DOUBLE_EQ(m.r2(10), 36.) << "r^2(10) = (0.05*10^2 + 1)^2 != 36.";
    ASSERT_DOUBLE_EQ(m.dr(5), 0.5) << "d/dx(0.05*x^2 + 1)|x=5 != 0.5";
    ASSERT_DOUBLE_EQ(m.dr(-5), 0.5) << "-d/dx(0.05*x^2 + 1)|x=-5 != 0.5";
}

TEST(ExponentialMesh, TestExponential1D)
{
    Exponential_mesh<1, double> m(1., 6., 0.02, 11);
    ASSERT_DOUBLE_EQ(m.r(10), 6.) << "r(10) = 0.02*(e^(0.02*10) - 1) + 1 != 6.";
    ASSERT_DOUBLE_EQ(m.r(-10), -4.) << "r(-10) = -0.02*(e^(0.02*10) - 1) + 1 != -4.";
    ASSERT_NEAR(m.r2(10), 36., TOL) << "r^2(10) = (0.02*(e^(0.02*10) - 1) + 1)^2 != 36.";
    ASSERT_NEAR(m.dr(5), 0.4991676378648055, TOL) << "d/dx(22.58..*(e^(0.02x) - 1) + 1)|x=5 != 0.4991676378648055";
    ASSERT_NEAR(m.dr(-5), 0.4991676378648055, TOL) << "-d/dx(22.58..*(e^(0.02x) - 1) + 1)|x=-5 != 0.4991676378648055";
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
            ASSERT_NEAR(dri, 0.4991676378648055, TOL) << "d/dx(22.58..*(e^(0.02x) - 1) + 1)|x=5 != 0.4991676378648055";
        });
    std::ranges::for_each(m.dr({-5, -5, -5, -5, -5}),
        [] (const double dri)
        {
            ASSERT_NEAR(dri, 0.4991676378648055, TOL) << "-d/dx(22.58..*(e^(0.02x) - 1) + 1)|x=-5 != 0.4991676378648055";
        });
}

TEST(MeshTest, TestIterate1D)
{
    double i = 0;
    Linear_mesh<1, double> m(0, 10, 11);
    for(auto v : m){
        ASSERT_NEAR(i++, v.r(), TOL);
    }
}

TEST(MeshTest, TestIterate5D)
{
    size_t i = 0;
    Linear_mesh<5, double> m({1., 1., 1., 1., 1.}, {6., 6., 6., 6., 6.}, {11, 11, 11, 11, 11});
    auto r_check{m.r()};
    for(auto ar : m){
        for(auto k = 0; k < 5; k++){
            ASSERT_NEAR(r_check[i][k], ar.r()[k], TOL) << "i " << i << ", k " << k;
        }
        i++;
    }
}
