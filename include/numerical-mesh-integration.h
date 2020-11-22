#ifndef NUMERICAL_MESH_LIB_INTEGRATION_H
#define NUMERICAL_MESH_LIB_INTEGRATION_H
#include "numerical-mesh.h"
#include <functional>
#include <type_traits>

template<class Scalar> requires std::is_floating_point<Scalar>::value
Scalar trapezoidal_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const Scalar& h, const std::function<Scalar(const Scalar&)>& integrand)
{
    Scalar res = 0;
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);
    Scalar distance = end - start;
    for(Scalar x = h; x <= distance; x += h){
        res +=  integrand(start + x - h)*mesh.dr(start + x - h) +
                integrand(start + x)*mesh.dr(start + x);
    }
    return res*h/2;
}

template<class Scalar> requires std::is_floating_point<Scalar>::value
Scalar simpson_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const Scalar& h, const std::function<Scalar(const Scalar&)>& integrand)
{
    Scalar res = 0;
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);
    Scalar distance = end - start;
    for(Scalar x = h; x <= distance/2; x += h){
        res +=  integrand(start + 2*x -2*h)*mesh.dr(start + 2*x - 2*h)
                 + 4*integrand(start + 2*x - h)*mesh.dr(start + 2*x - h)
                 + integrand(start + 2*x)*mesh.dr(start + 2*x);
    }
    return res*h/3;
}

template<class Scalar, size_t K> requires std::is_floating_point<Scalar>::value
std::array<Scalar, K> end_point_corrections()
{
    return {9./24, 28./24, 23./24};
}

template<class Scalar, size_t K = 3> requires std::is_floating_point<Scalar>::value
Scalar corrected_trapezoidal_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const Scalar& h, const std::function<Scalar(const Scalar&)>& integrand)
{
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);
    Scalar distance = end - start;
    if(distance < 2*K*h){
        return 0;
    }
    auto weights = end_point_corrections<Scalar, K>();
    Scalar res = 0;
    size_t i = 0;
    for(Scalar x = 0; x < h*K; x += h){
        res += weights[i]*integrand(start + x)*mesh.dr(start + x);
        res += weights[i]*integrand(end - x)*mesh.dr(end - x);
        i++;
    }
    for(Scalar x = K*h; x <= distance - K*h; x += h){
        res += integrand(start + x)*mesh.dr(start + x);
    }
    return res*h;
}



#endif  //NUMERICAL_MESH_LIB_INTEGRATION_H
