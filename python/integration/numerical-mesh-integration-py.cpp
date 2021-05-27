#include <numerical-mesh-integration.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(integration, m) {
    m.attr("__name__") = "NumericalMesh.Integration";

    m.def("trapezoid", static_cast<double (*)(const Mesh_base<1, double>&, const std::function<double(const double&)>&)>
        (&corrected_trapezoidal_integral<double, 3>));
    m.def("correctedTrapezoid", static_cast<double (*)(const Mesh_base<1, double>&, const std::function<double(const double&)>&)>
        (&corrected_trapezoidal_integral<double, 3>));
    m.def("simpson", static_cast<double (*)(const Mesh_base<1, double>&, const std::function<double(const double&)>&)>
        (&corrected_trapezoidal_integral<double, 3>));
}
