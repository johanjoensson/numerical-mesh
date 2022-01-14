#include <numerical-mesh-integration.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>


namespace py = pybind11;

void init_integration(py::module& m){
    py::module_ m_integration = m.def_submodule("Integration", "Numerical integration on NumericalMesh");
    m_integration.def("trapezoid", static_cast<double (*)(const Mesh_base<1, double>&, const std::function<double(const Mesh_base<1, double>::mesh_point&)>&)>
        (&trapezoidal_integral<double>));
    m_integration.def("correctedTrapezoid", static_cast<double (*)(const Mesh_base<1, double>&, const std::function<double(const Mesh_base<1, double>::mesh_point&)>&)>
        (&corrected_trapezoidal_integral<double, 3>));
    m_integration.def("simpson", static_cast<double (*)(const Mesh_base<1, double>&, const std::function<double(const Mesh_base<1, double>::mesh_point&)>&)>
        (&simpson_integral<double>));
}
