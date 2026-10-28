/**
 * @file kesmodel_py.cpp
 * @brief Python module kesmodel: the C++ key equation solvers with the interface of kes_ref.py.
 *
 * Snapshots and traces are returned as dicts with the same keys as kes_ref, so traces of both
 * models can be compared with ==.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "gf.hpp"
#include "kes_core.hpp"

namespace py = pybind11;
using fec::GaloisField;
using fec::Kes;
using fec::KesResult;
using fec::Snapshot;
using fec::Symbol;

namespace {

py::dict to_dict(const Snapshot& s) {
    py::dict d;
    d["r"] = s.r;
    d["gamma"] = s.gamma;
    d["k"] = s.k;
    d["mc"] = s.mc;
    d["discrepancy"] = s.discrepancy;
    for (const auto& [name, values] : s.registers) {
        d[py::str(name)] = py::cast(values);
    }
    return d;
}

[[noreturn]] void zero_division(const char* message) {
    PyErr_SetString(PyExc_ZeroDivisionError, message);
    throw py::error_already_set();
}

}  // namespace

PYBIND11_MODULE(kesmodel, m) {
    m.doc() = "C++ register-level models of the key equation solvers (see kes_ref.py)";

    py::class_<GaloisField>(m, "GaloisField")
        .def(py::init<unsigned, std::uint32_t>(), py::arg("m"), py::arg("poly"))
        .def_property_readonly("m", &GaloisField::m)
        .def_property_readonly("poly", &GaloisField::poly)
        .def_property_readonly("order", &GaloisField::order)
        .def_property_readonly("n", &GaloisField::n)
        .def("mul", [](const GaloisField& f, Symbol a, Symbol b) { return f.mul(f.check(a), f.check(b)); })
        .def("inverse", [](const GaloisField& f, Symbol a) {
            if (f.check(a) == 0) {
                zero_division("0 has no multiplicative inverse");
            }
            return f.inverse(a);
        })
        .def("exp", &GaloisField::exp)
        .def("log", [](const GaloisField& f, Symbol a) {
            if (f.check(a) == 0) {
                zero_division("log(0) is undefined");
            }
            return f.log(a);
        })
        .def("__repr__", [](const GaloisField& f) {
            return "kesmodel.GaloisField(m=" + std::to_string(f.m()) + ", poly=" + std::to_string(f.poly()) + ")";
        });

    py::class_<KesResult>(m, "KesResult")
        .def_readonly("lam", &KesResult::lam)
        .def_readonly("omega", &KesResult::omega)
        .def_readonly("form", &KesResult::form)
        .def_readonly("iterations", &KesResult::iterations)
        .def_readonly("length", &KesResult::length)
        .def_property_readonly("degree", &KesResult::degree)
        .def_property_readonly("trace", [](const KesResult& res) -> py::object {
            if (!res.traced) {
                return py::none();
            }
            py::list states;
            for (const auto& s : res.trace) {
                states.append(to_dict(s));
            }
            return states;
        });

    py::class_<Kes>(m, "Kes")
        .def_property_readonly("name", &Kes::name)
        .def_property_readonly("form", &Kes::form)
        .def_property_readonly("iterations", &Kes::iterations)
        .def_property_readonly("t", &Kes::t)
        .def_property_readonly("r", &Kes::r)
        .def_property_readonly("loaded", &Kes::loaded)
        .def_property_readonly("done", &Kes::done)
        .def_property_readonly("gamma", &Kes::gamma)
        .def_property_readonly("k", &Kes::k)
        .def_property_readonly("mc", &Kes::mc)
        .def_property_readonly("discrepancy", &Kes::discrepancy)
        .def("reset", &Kes::reset)
        .def("load", &Kes::load, py::arg("syndromes"))
        .def("step", &Kes::step)
        .def("result", &Kes::result)
        .def("run", &Kes::run, py::arg("syndromes"), py::arg("trace") = false)
        .def("snapshot", [](const Kes& kes) { return to_dict(kes.snapshot()); });

    m.def("make_kes", &fec::make_kes, py::arg("name"), py::arg("field"), py::arg("t"),
          "Build a solver: 'ibm', 'ribm' or 'ribm_systolic'.");
    m.attr("ARCHITECTURES") = py::make_tuple("ibm", "ribm", "ribm_systolic");
}
