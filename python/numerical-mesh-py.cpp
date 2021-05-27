#include <numerical-mesh.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

template<size_t dim, arithmetic Scalar>
class PyMesh : public Mesh_base<dim, Scalar> {
public:

    using ClassType = Mesh_base<dim, Scalar>;
    using Mesh_base<dim, Scalar>::Mesh_base;

    Scalar r(const Scalar& x) const noexcept override
    {

        PYBIND11_OVERRIDE_PURE(
             Scalar,
             ClassType,
             r,
             x
        );
    }

    Scalar r2(const Scalar& x) const noexcept override
    {
        PYBIND11_OVERRIDE_PURE(
            Scalar,                 /* Return type */
            ClassType,   /* Parent class */
            r2,                      /* Name of function in C++ (must match Python name) */
            x                       /* Argument(s) */
        );
    }

    Scalar dr(const Scalar& x) const noexcept override
    {
        PYBIND11_OVERRIDE_PURE(
            Scalar,                 /* Return type */
            ClassType,   /* Parent class */
            dr,                      /* Name of function in C++ (must match Python name) */
            x                       /* Argument(s) */
        );
    }
    Scalar d2r(const Scalar& x) const noexcept override
    {
        PYBIND11_OVERRIDE_PURE(
            Scalar,                 /* Return type */
            ClassType,   /* Parent class */
            d2r,                      /* Name of function in C++ (must match Python name) */
            x                       /* Argument(s) */
        );
    }
    Scalar d3r(const Scalar& x) const noexcept override
    {
        PYBIND11_OVERRIDE_PURE(
            Scalar,                 /* Return type */
            ClassType,   /* Parent class */
            d3r,                      /* Name of function in C++ (must match Python name) */
            x                       /* Argument(s) */
        );
    }

};

PYBIND11_MODULE(NumericalMesh, m) {
    m.doc() = "pybind11 Python bindings for the numerical mesh C++ library.";

    py::class_<Mesh_base<1, double>, PyMesh<1, double>> (m, "MeshBase")
        .def(py::init<double, double, double, size_t>())
        .def("r", static_cast<std::vector<double> (Mesh_base<1, double>::*)()const> (&Mesh_base<1, double>::r))
        .def("r2", static_cast<std::vector<double> (Mesh_base<1, double>::*)()const>(&Mesh_base<1, double>::r2))
        .def("dr", static_cast<std::vector<double> (Mesh_base<1, double>::*)()const>(&Mesh_base<1, double>::dr))
        .def("d2r", static_cast<std::vector<double> (Mesh_base<1, double>::*)()const>(&Mesh_base<1, double>::d2r))
        .def("d3r", static_cast<std::vector<double> (Mesh_base<1, double>::*)()const>(&Mesh_base<1, double>::d3r))
        .def("r", static_cast<double (Mesh_base<1, double>::*)(const double&)const> (&Mesh_base<1, double>::r))
        .def("r2", static_cast<double (Mesh_base<1, double>::*)(const double&)const>(&Mesh_base<1, double>::r2))
        .def("dr", static_cast<double (Mesh_base<1, double>::*)(const double&)const>(&Mesh_base<1, double>::dr))
        .def("d2r", static_cast<double (Mesh_base<1, double>::*)(const double&)const>(&Mesh_base<1, double>::d2r))
        .def("d3r", static_cast<double (Mesh_base<1, double>::*)(const double&)const>(&Mesh_base<1, double>::d3r))
        .def("size", &Mesh_base<1, double>::size)
        .def("dim", &Mesh_base<1, double>::dim)
        .def("__len__", &Mesh_base<1, double>::size)
        .def("__iter__", [](const Mesh_base<1, double> &m)
            { return py::make_iterator(m.begin(), m.end()); }, py::keep_alive<0, 1>())
    ;

    py::class_<Mesh_base<1, double>::mesh_point> (m, "MeshPoint")
        .def("r", &Mesh_base<1, double>::mesh_point::r)
        .def("r2", &Mesh_base<1, double>::mesh_point::r2)
        .def("dr", &Mesh_base<1, double>::mesh_point::dr)
        .def("d2r", &Mesh_base<1, double>::mesh_point::d2r)
        .def("d3r", &Mesh_base<1, double>::mesh_point::d3r)
        ;

    py::class_<Linear_mesh<1, double>, Mesh_base<1, double>>(m, "LinearMesh")
        .def(py::init<double, double, size_t>())
    ;
    py::class_<Quadratic_mesh<1, double>, Mesh_base<1, double>>(m, "QuadraticMesh")
        .def(py::init<double, double, size_t>())
    ;
    py::class_<Exponential_mesh<1, double>, Mesh_base<1, double>>(m, "ExponentialMesh")
        .def(py::init<double, double, double, size_t>())
    ;
}
