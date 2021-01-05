#ifndef NUMERICAL_MESH_LIB_H
#define NUMERICAL_MESH_LIB_H

#if __cplusplus < 201103L // Require at least c++11
#error The numerical-mesh library requires at least c++11!
#endif

#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <iostream>

/*******************************************************************************
* The __cplusplus macro is used to determine which c++ standard is used for the*
* compilation. This is "needed" because I want to use this library in projects *
* using c++11 as well as projects using c++20. This means parts of the code had*
* to be rewritten in order to enable the cool c++20 features only if c++20 is  *
* used. As a side effect this also made it possible to write code for c++17 as *
* well.                                                                        *
* Features requiring c++17:                                                    *
*  -- structured bindings i.e., auto [a,b,c,...] = ...                         *
* Features requiring c++20:                                                    *
*  -- Concepts e.g., integral and floating_point                               *
*  -- ranges-versions of algorithms                                            *
*  -- std::transform with 2 input containers                                   *
*******************************************************************************/

/*******************************************************************************
* The arithmetic concept, integral or floating point, is used to signify any
* numeric type.
*******************************************************************************/
#if __cplusplus >= 202002L
concept arithmetic = std::integral || std::floating_point
#endif

namespace{
#if __cplusplus >= 202002L
    template<arithmetic T>
#else
    template<typename T>
#endif
    constexpr int sgn(T val)
    {
    return (T(0) < val) - (val < T(0));
    }
}

/*******************************************************************************
* Base class for all numerical mesh classes. All meshes can be described by
* (max) three parameters and the number of mesh points (including endpoints).
* @param alpha First parameter for the mesh.
* @param beta Second parameter for the mesh.
* @param gamma Third parameter for the mesh.
* @param N Number of mesh points in the mesh (including endpoints).
*******************************************************************************/
#if __cplusplus >= 202002L
template<size_t Dim, floating_point Scalar = double>
#else
template<size_t Dim, class Scalar = double>
#endif
class Mesh_base {
    using Arr = std::array<Scalar, Dim>;

protected:
    template<class A, class B>
    std::array<std::tuple<A, B>, Dim> zip2 (const std::array<A, Dim>& a, const std::array<B, Dim>& b) const
    {
        std::array<std::tuple<A, B>, Dim> res;
#if __cplusplus >= 202002L
        std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(res),
            [](const A& ai, const B& bi)
            {
                return std::tuple<A, B>(ai, bi);
            });
#else
        auto ai = a.begin();
        auto bi = b.begin();
        auto ri = res.begin();
        for(;ai != a.end(); ai++, bi++, ri++){
            *ri = {*ai, *bi};
        }
#endif
            return res;
    }

    template<class A, class B, class C>
    std::array<std::tuple<A, B, C>, Dim> zip3(const std::array<A, Dim>& a, const std::array<B, Dim>& b, std::array<C, Dim> c) const
    {
        auto tmp = zip2(b, c);
        std::array<std::tuple<A, B, C>, Dim> res;
#if __cplusplus >= 202002L
        std::transform(std::begin(a), std::end(a), std::begin(tmp), std::begin(res),
            [] (const Scalar ai, const std::tuple<B, C> ti)
            {
                auto [bi, ci] = ti;
                return std::tuple<A, B, C>(ai, bi, ci);
            });
#else
        auto ai = a.begin();
        auto bi = b.begin();
        auto ci = c.begin();
        auto ri = res.begin();
        for(;ai != a.end(); ai++, bi++, ci++, ri++){
            *ri = {*ai, *bi, *ci};
        }
#endif
        return res;
    }

    template<class A, class B, class C, class D>
    std::array<std::tuple<A, B, C, D>, Dim> zip4 (const std::array<A, Dim>& a, const std::array<B, Dim>& b, const std::array<C, Dim>& c, const std::array<D, Dim>& d) const
    {
        auto z1 = zip2(a, b), z2 = zip2(c, d);
        std::array<std::tuple<A, B, C, D>, Dim> res;
#if __cplusplus >= 202002L
        std::transform(std::begin(z1), std::end(z1), std::begin(z2), std::begin(res),
        [] (const std::tuple<A, B> ti, const std::tuple<C, D> ui)
        {
            auto [ai, bi] = ti;
            auto [ci, di] = ui;
            return std::tuple<A, B, C, D>(ai, bi, ci, di);
        });
#else
        auto ai = a.begin();
        auto bi = b.begin();
        auto ci = c.begin();
        auto di = d.begin();
        auto ri = res.begin();
        for(;ai != a.end(); ai++, bi++, ci++, di++, ri++){
            *ri = {*ai, *bi, *ci, *di};
        }
#endif
        return res;
    }

    template<class A, class B, class C, class D>
    std::array<std::tuple<A, B, C, D>, Dim> zip4 (const std::array<A, Dim>& a, const std::array<std::tuple<B, C, D>, Dim>& t) const
    {
        std::array<std::tuple<A, B, C, D>, Dim> res;
#if __cplusplus >= 202002L
        std::transform(std::begin(a), std::end(a), std::begin(t), std::begin(res),
        [] (const A ai, const std::tuple<B, C, D> ti)
        {
            auto [bi, ci, di] = ti;
            return std::tuple<A, B, C, D>(ai, bi, ci, di);
        });
#else
        auto ai = a.begin();
        auto ti = t.begin();
        auto ri = res.begin();
        for(;ai != a.end(); ai++, ti++, ri++){
            auto bi = std::get<0>(*ti);
            auto ci = std::get<1>(*ti);
            auto di = std::get<2>(*ti);
            *ri = {*ai, bi, ci, di};
        }
#endif
        return res;
    }

    Arr alpha_m, beta_m, gamma_m;
    std::array<size_t, Dim> N_m;
    Mesh_base() = default;

public:
    virtual ~Mesh_base() = default;
    Mesh_base(const Arr alpha, const Arr beta, const Arr gamma, const std::array<size_t, Dim>& N)
     : alpha_m(alpha), beta_m(beta), gamma_m(gamma), N_m(N)
    {}

    Mesh_base(const Mesh_base&) = default;
    Mesh_base(Mesh_base&&) = default;

    Mesh_base& operator=(const Mesh_base&) = default;
    Mesh_base& operator=(Mesh_base&&) = default;

    auto parameters() const {return zip3(alpha_m, beta_m, gamma_m);}

    virtual Arr r(const Arr& x) const = 0;
    virtual Arr r2(const Arr& x) const = 0;
    virtual Arr dr(const Arr& x) const = 0;

    Arr operator()(const Arr& x) const {return r(x);}

    std::array<std::vector<Scalar>, Dim> r() const;

    std::array<size_t, Dim> dim() const {return N_m;}

};

/*******************************************************************************
* 1D template specialization, instead of using 1D-arrays we can use the vales
* directly.
* @param alpha First parameter for the mesh.
* @param beta Second parameter for the mesh.
* @param gamma Third parameter for the mesh.
* @param N Number of mesh points in the mesh (including endpoints).
*******************************************************************************/
#if __cplusplus >= 202002L
template<floating_point Scalar>
#else
template<class Scalar>
#endif
class Mesh_base<1, Scalar> {
protected:
    Scalar alpha_m, beta_m, gamma_m;
    size_t N_m;
    Mesh_base() = default;

public:
    virtual ~Mesh_base() = default;
    Mesh_base(const Scalar alpha, const Scalar beta, const Scalar gamma, const size_t N)
     : alpha_m(alpha), beta_m(beta), gamma_m(gamma), N_m(N)
    {}

    Mesh_base(const Mesh_base&) = default;
    Mesh_base(Mesh_base&&) = default;

    Mesh_base& operator=(const Mesh_base&) = default;
    Mesh_base& operator=(Mesh_base&&) = default;

    virtual Scalar r(const Scalar& x) const = 0;
    virtual Scalar r2(const Scalar& x) const = 0;
    virtual Scalar dr(const Scalar& x) const = 0;
    Scalar operator()(const Scalar& x) const {return r(x);}

    virtual std::vector<Scalar> r() const
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->r(i));
        }
        return res;
    }

    virtual std::vector<Scalar> r2() const
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->r2(i));
        }
        return res;
    }

    virtual std::vector<Scalar> dr() const
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->dr(i));
        }
        return res;
    }

    template<class Func>
    std::vector<Scalar> evaluate(Func f, size_t h = 1)
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i += h){
            res.push_back(f(this->r(i)));
        }
        return res;
    }

    template<class Func>
    std::vector<Scalar> evaluate(Func f, Scalar x0, Scalar x1, Scalar h)
    {
        std::vector<Scalar> res;
        for(Scalar x = x0; x <= x1; x += h){
            res.push_back(f(this->r(x)));
        }
        return res;
    }

    size_t dim() const {return N_m;}

    class iterator{
    private:
        const Mesh_base& m_m;
        Scalar i_m;
    public:
        typedef std::random_access_iterator_tag iterator_category;
        iterator() = delete;
        iterator(const iterator&) = default;
        iterator(iterator&&) = default;
        iterator(const Mesh_base& m, const Scalar i) : m_m(m), i_m(i){}
        ~iterator() = default;
        iterator& operator=(const iterator&) = default;
        iterator& operator=(iterator&&) = default;

        Scalar operator*(){return m_m.r(i_m);}

        bool operator==(const iterator& b)const{return i_m == b.i_m;}
        bool operator!=(const iterator& b)const{return !(*this == b);}
        bool operator>=(const iterator& b)const{return i_m >= b.i_m;}
        bool operator<=(const iterator& b)const{return i_m <= b.i_m;}
        bool operator>(const iterator& b)const{return !(*this <= b);}
        bool operator<(const iterator& b)const{return !(*this >= b);}

        iterator operator++(int){auto res = *this; i_m++; return res;}
        iterator& operator++(){i_m++; return *this;}
        iterator operator--(int){auto res = *this; i_m--; return res;}
        iterator& operator--(){i_m--; return *this;}

        iterator& operator+=(const size_t n){i_m += n; return *this;}
        iterator& operator-=(const size_t n){i_m -= n; return *this;}

        iterator operator+(const size_t n)const{auto res = *this; return (res += n);}
        iterator operator-(const size_t n)const{auto res = *this; return (res -= n);}
        size_t operator-(const iterator& it)const{return this->i_m - it.i_m;}

        Scalar operator[](const size_t n)const{return m_m.r(n);}
    };
    class const_iterator{
    private:
        const Mesh_base& m_m;
        Scalar i_m;
    public:
        typedef std::random_access_iterator_tag iterator_category;
        const_iterator() = delete;
        const_iterator(const const_iterator&) = default;
        const_iterator(const_iterator&&) = default;
        const_iterator(const Mesh_base& m, const Scalar i) : m_m(m), i_m(i){}
        ~const_iterator() = default;
        const_iterator& operator=(const const_iterator&) = default;
        const_iterator& operator=(const_iterator&&) = default;

        Scalar operator*(){return m_m.r(i_m);}

        bool operator==(const const_iterator& b)const{return i_m == b.i_m;}
        bool operator!=(const const_iterator& b)const{return !(*this == b);}
        bool operator>=(const const_iterator& b)const{return i_m >= b.i_m;}
        bool operator<=(const const_iterator& b)const{return i_m <= b.i_m;}
        bool operator>(const const_iterator& b)const{return !(*this <= b);}
        bool operator<(const const_iterator& b)const{return !(*this >= b);}

        const_iterator operator++(int){auto res = *this; i_m++; return res;}
        const_iterator& operator++(){i_m++; return *this;}
        const_iterator operator--(int){auto res = *this; i_m--; return res;}
        const_iterator& operator--(){i_m--; return *this;}

        const_iterator& operator+=(const size_t n){i_m += n; return *this;}
        const_iterator& operator-=(const size_t n){i_m -= n; return *this;}

        const_iterator operator+(const size_t n)const{auto res = *this; return (res += n);}
        const_iterator operator-(const size_t n)const{auto res = *this; return (res -= n);}
        size_t operator-(const const_iterator& it)const{return this->i_m - it.i_m;}

        const Scalar operator[](const size_t n)const{return m_m.r(n);}
    };

    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    iterator begin(){return iterator(*this, 0);}
    const_iterator begin() const {return const_iterator(*this, 0);}
    const_iterator cbegin() const {return const_iterator(*this, 0);}
    iterator end(){return iterator(*this, N_m);}
    const_iterator end() const {return const_iterator(*this, N_m);}
    const_iterator cend() const {return const_iterator(*this, N_m);}

    reverse_iterator rbegin(){return reverse_iterator(this->begin());}
    const_reverse_iterator rbegin() const {return const_reverse_iterator(this->begin());}
    const_reverse_iterator crbegin() const {return const_reverse_iterator(this->begin());}
    reverse_iterator rend(){return reverse_iterator(this->end());}
    const_reverse_iterator rend() const {return const_reverse_iterator(this->end());}
    const_reverse_iterator crend() const {return const_reverse_iterator(this->cend());}
};

/*******************************************************************************
* Linear mesh in arbitrary dimensions.
* @param R_min Array containing the lower end points for the mesh.
* @param R_max Array containing the upper end points for the mesh.
* @param N Array containing the number of mesh points along the different
*          dimensions of the mesh.
* @param alpha Array containing the step sizes along the different dimensions.
* @param beta Not used in a linear mesh.
* @param gamma R_min
*******************************************************************************/
template<size_t D, class Scalar = double>
class Linear_mesh : public Mesh_base<D, Scalar> {
private:
    Linear_mesh() = default;

public:
    using Arr = std::array<Scalar, D>;
    Linear_mesh(const Arr& R_min, const Arr& R_max, const std::array<size_t, D>& N)
    : Mesh_base<D, Scalar>(
        [this] (const Arr& R_min, const Arr& R_max,  const std::array<size_t, D>& N)
        {
            Arr res;
            auto z3 = this->zip3(R_min,  R_max, N);
            std::transform(std::begin(z3), std::end(z3), std::begin(res),
                [] (std::tuple<Scalar, Scalar, size_t> ti)
                {
#if __cplusplus >= 201703L
                    auto [r1, r2, n] = ti;
#else
                    auto r1 = std::get<0>(ti);
                    auto r2 = std::get<1>(ti);
                    auto n = std::get<2>(ti);
#endif
                    return (r2 - r1)/(n - 1);
                });
                return res;
            }(R_min, R_max, N)
        , Arr(), R_min, N)
    {}

    ~Linear_mesh() = default;

    Linear_mesh(const Linear_mesh&) = default;
    Linear_mesh(Linear_mesh&&) = default;

    Linear_mesh& operator=(const Linear_mesh&) = default;
    Linear_mesh& operator=(Linear_mesh&&) = default;

    Arr r(const Arr& x) const override
    {
        Arr res;

#if __cplusplus >= 202002L
        std::transform(std::begin(x), std::end(x), std::begin(this->parameters()), std::begin(res),
        [] (const Scalar xi, const std::tuple<Scalar, Scalar, Scalar> ti)
        {
            auto [ai, bi, ci] = ti;
            return ai*xi +0*bi + ci;
        });
#else
        auto z4 = this->zip4(x, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
        [] (const std::tuple<Scalar, Scalar, Scalar, Scalar> ti)
        {
#if __cplusplus >= 201703L
            auto [xi, ai, bi, ci] = ti;
#else
            auto xi = std::get<0>(ti);
            auto ai = std::get<1>(ti);
            auto bi = std::get<2>(ti);
            auto ci = std::get<3>(ti);
#endif
            return ai*xi +0*bi + ci;
        });
#endif
        return res;
    }

    Arr r2(const Arr& x) const override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res),
            [] (const Scalar s) {return s*s;});
        return res;
    }

    Arr dr(const Arr& x) const override
    {
        return this->alpha_m;
    }
};

/*******************************************************************************
* Linear mesh in one dimension, r = alpha*x + gamma.
* @param R_min Lower end point for the mesh.
* @param R_max Upper end point for the mesh.
* @param N Number of mesh points along the mesh.
* @param alpha Step size.
* @param beta Not used in a linear mesh.
* @param gamma R_min
*******************************************************************************/
template<class Scalar>
class Linear_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
private:
    Linear_mesh() = default;

public:
    Linear_mesh(const Scalar R_min, const Scalar R_max, const size_t N)
    : Mesh_base<1, Scalar>((R_max - R_min)/(N - 1), 0, R_min, N)
    {}

    ~Linear_mesh() = default;

    Linear_mesh(const Linear_mesh&) = default;
    Linear_mesh(Linear_mesh&&) = default;

    Linear_mesh& operator=(const Linear_mesh&) = default;
    Linear_mesh& operator=(Linear_mesh&&) = default;

    using Mesh_base<1, Scalar>::r;
    using Mesh_base<1, Scalar>::r2;
    using Mesh_base<1, Scalar>::dr;

    Scalar r(const Scalar& x) const override
    {
        return x*this->alpha_m + this->gamma_m;
    }

    Scalar r2(const Scalar& x) const override
    {
        return std::pow(this->r(x), 2);
    }

    Scalar dr(const Scalar& x) const override
    {
        return this->alpha_m;
    }
};

/*******************************************************************************
* Quadratic mesh in arbitrary dimensions, r_i = alpha_i*x_i*x_i + gamma_i.
* @param R_min Array containing the lower end points for the mesh.
* @param R_max Array containing the upper end points for the mesh.
* @param N Array containing the number of mesh points along the different
*          dimensions of the mesh.
* @param alpha Array containing the first mesh parameter for each dimension.
* @param beta Not used in a quadratic mesh.
* @param gamma R_min
*******************************************************************************/
template<size_t D, class Scalar = double>
class Quadratic_mesh : public Mesh_base<D, Scalar> {
public:
    using Arr = std::array<Scalar, D>;
    Quadratic_mesh() = delete;
    Quadratic_mesh(const Arr& R_min, const Arr& R_max, const std::array<size_t, D>& N)
    : Mesh_base<D, Scalar>(
        [this] (const Arr& R_min, const Arr& R_max, const std::array<size_t, D>& N)
        {
            Arr res;
            auto z3 = this->zip3(R_min, R_max, N);
            std::transform(std::begin(z3), std::end(z3), std::begin(res),
                [] (const std::tuple<Scalar, Scalar, size_t>& ti)
                {
#if __cplusplus >= 201703L
                    auto [r1, r2, n] = ti;
#else
                    auto r1 = std::get<0>(ti);
                    auto r2 = std::get<1>(ti);
                    auto n = std::get<2>(ti);
#endif
                    return (r2-r1)/std::pow(n - 1, 2);
                });
            return res;
        }(R_min, R_max, N)
        , Arr(), R_min, N )
    {}

    ~Quadratic_mesh() = default;

    Quadratic_mesh(const Quadratic_mesh&) = default;
    Quadratic_mesh(Quadratic_mesh&&) = default;

    Quadratic_mesh& operator=(const Quadratic_mesh&) = default;
    Quadratic_mesh& operator=(Quadratic_mesh&&) = default;

    Arr r(const Arr& x) const override
    {
        Arr res;
#if __cplusplus >= 202002L
        std::transform(std::begin(x), std::end(x), std::begin(this->parameters()), std::begin(res),
            [] (const Scalar& xi, const std::tuple<Scalar, Scalar, Scalar>& ti)
            {
                auto [a, b, c] = ti;
                return sgn(xi)*a*std::pow(xi, 2) + c;
            });
#else
        auto z4 = this->zip4(x, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
            [] (const std::tuple<Scalar, Scalar, Scalar, Scalar>& ti)
            {
#if __cplusplus >= 201703L
                auto [xi, a, b, c] = ti;
#else
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
#endif
                return sgn(xi)*a*std::pow(xi, 2) + c;
            });
#endif
        return res;
    }

    Arr r2(const Arr& x) const override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res), [] (const Scalar& s) {return s*s;});
        return res;
    }

    Arr dr(const Arr& x) const override
    {
        Arr res;
#if __cplusplus >= 202002L
        std::transform(std::begin(x), std::end(x), std::begin(this->parameters()), std::begin(res),
            [] (const Scalar xi, const std::tuple<Scalar, Scalar, Scalar> ti)
            {
                auto [a, b, c] = ti;
                return sgn(xi)*2*a*xi + b*0 + c*0;
            });
#else
        auto z4 = this->zip4(x, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
            [] (const std::tuple<Scalar, Scalar, Scalar, Scalar>& ti)
            {
#if __cplusplus >= 201703L
                auto [xi, a, b, c] = ti;
#else
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
#endif
                return sgn(xi)*2*a*xi + b*0 + c*0;
            });
#endif
        return res;
    }
};

/*******************************************************************************
* Quadratic mesh in one dimension, r = alpha*x*x + gamma.
* @param R_min Lower end point for the mesh.
* @param R_max Upper end point for the mesh.
* @param N Number of mesh points along the mesh.
* @param alpha Array containing the first mesh parameter.
* @param beta Not used in a quadratic mesh.
* @param gamma R_min
*******************************************************************************/
template<class Scalar>
class Quadratic_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
public:
    Quadratic_mesh() = delete;
    Quadratic_mesh(const Scalar R_min, const Scalar R_max, const size_t N)
    : Mesh_base<1, Scalar>((R_max - R_min)/std::pow(N - 1, 2), 0, R_min, N )
    {}

    ~Quadratic_mesh() = default;

    Quadratic_mesh(const Quadratic_mesh&) = default;
    Quadratic_mesh(Quadratic_mesh&&) = default;

    Quadratic_mesh& operator=(const Quadratic_mesh&) = default;
    Quadratic_mesh& operator=(Quadratic_mesh&&) = default;

    Scalar r(const Scalar& x) const override
    {
        return sgn(x)*std::pow(x, 2)*this->alpha_m + this->gamma_m;
    }

    Scalar r2(const Scalar& x) const override
    {
        return std::pow(this->r(x), 2);
    }

    Scalar dr(const Scalar& x) const override
    {
        return 2*this->alpha_m*sgn(x)*x;
    }
};

/*******************************************************************************
* Exponential mesh in arbitrary dimensions, r_i = alpha_i(exp(beta_i*x) - 1) +
* gamma_i.
* @param R_min Array containing the lower end points for the mesh.
* @param R_max Array containing the upper end points for the mesh.
* @param N Array containing the number of mesh points along the dimensions of
*          the mesh.
* @param alpha Array containing the first mesh parameters.
* @param beta Array containing the second mesh parameters.
* @param gamma R_min.
*******************************************************************************/
template<size_t D, class Scalar = double>
class Exponential_mesh : public Mesh_base<D, Scalar> {
public:
    using Arr = std::array<Scalar, D>;
    Exponential_mesh() = delete;
    Exponential_mesh(const Arr& R_min, const Arr& R_max, const Arr& beta, const std::array<size_t, D>& N)
    : Mesh_base<D, Scalar>([this] (const Arr& R_min, const Arr& R_max, const Arr& beta, const std::array<size_t, D>& N)
        {
            Arr res;
            auto z4 = this->zip4(R_min, R_max, beta, N);
            std::transform(std::begin(z4), std::end(z4), std::begin(res),
                [] (const std::tuple<Scalar, Scalar, Scalar, size_t> ti)
                {
#if __cplusplus >= 201703L
                    auto [r1, r2, b, n] = ti;
#else
                    auto r1 = std::get<0>(ti);
                    auto r2 = std::get<1>(ti);
                    auto b = std::get<2>(ti);
                    auto n = std::get<2>(ti);
#endif
                    return (r2 - r1)/(std::exp(b*(n - 1)) - 1);
                });
            return res;
        }(R_min, R_max, beta, N), beta, R_min, N)
    {}

    ~Exponential_mesh() = default;

    Exponential_mesh(const Exponential_mesh&) = default;
    Exponential_mesh(Exponential_mesh&&) = default;

    Exponential_mesh& operator=(const Exponential_mesh&) = default;
    Exponential_mesh& operator=(Exponential_mesh&&) = default;

    Arr r(const Arr& x) const override
    {
        Arr res;
#if __cplusplus >= 202002L
        std::transform(std::begin(x), std::end(x), std::begin(this->parameters()), std::begin(res),
            [] (const Scalar& xi, const std::tuple<Scalar, Scalar, Scalar>& ti)
            {
                auto [a, b, c] = ti;
                return sgn(xi)*a*(std::exp(sgn(xi)*xi*b) - 1) + c;
            });
#else
        auto z4 = this->zip4(x, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
            [] (const std::tuple<Scalar, Scalar, Scalar, Scalar>& ti)
            {
#if __cplusplus >= 201703L
                auto [xi, a, b, c] = ti;
#else
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
#endif
                return sgn(xi)*a*(std::exp(sgn(xi)*xi*b) - 1) + c;
            });
#endif
        return res;
    }

    Arr r2(const Arr& x) const override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res), [] (const Scalar& s) {return s*s;});
        return res;
    }

    Arr dr(const Arr& x) const override
    {
        Arr x_abs;
        std::transform(std::begin(x), std::end(x), std::begin(x_abs),
            [] (const Scalar& xi)
            {
                return sgn(xi)*xi;
            }
        );
        Arr r = this->r(x_abs);

        Arr res;
#if __cplusplus >= 202002L
        std::transform(std::begin(r), std::end(r), std::begin(this->parameters()), std::begin(res),
            [] (const Scalar& ri, const std::tuple<Scalar, Scalar, Scalar>& ti)
            {
                auto [a, b, c] = ti;
                return b*(ri + a - c);
            });
#else
        auto z4 = this->zip4(r, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
            [] (const std::tuple<Scalar, Scalar, Scalar, Scalar>& ti)
            {
#if __cplusplus >= 201703L
                auto [ri, a, b, c] = ti;
#else
                auto ri = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
#endif
                return b*(ri + a - c);
            });
#endif
        return res;
    }
};

/*******************************************************************************
* Exponential mesh in one dimension, r = alpha(exp(beta*x) - 1) + gamma.
* @param R_min Lower end point for the mesh.
* @param R_max Upper end point for the mesh.
* @param N Number of mesh points along the dimensions of the mesh.
* @param alpha First mesh parameter.
* @param beta Second mesh parameter.
* @param gamma R_min.
*******************************************************************************/
template<class Scalar>
class Exponential_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
public:
    Exponential_mesh() = delete;
    Exponential_mesh(const Scalar R_min, const Scalar R_max, const Scalar beta, const size_t N)
    : Mesh_base<1, Scalar>((R_max - R_min)/std::expm1(beta*(N-1)), beta, R_min, N)
    {}

    ~Exponential_mesh() = default;

    Exponential_mesh(const Exponential_mesh&) = default;
    Exponential_mesh(Exponential_mesh&&) = default;

    Exponential_mesh& operator=(const Exponential_mesh&) = default;
    Exponential_mesh& operator=(Exponential_mesh&&) = default;

    Scalar r(const Scalar& x) const override
    {
        Scalar a = this->alpha_m, b = this->beta_m, c = this->gamma_m;
        return sgn(x)*a*std::expm1(sgn(x)*b*x) + c;
    }

    Scalar r2(const Scalar& x) const override
    {
        return std::pow(this->r(x), 2);
    }

    Scalar dr(const Scalar& x) const override
    {
        Scalar a = this->alpha_m, b = this->beta_m, c = this->gamma_m;
        return b*(this->r(sgn(x)*x) + a - c);
    }
};

#endif //NUMERICAL_MESH_LIB_H
