/**
 * @file kessc_py.cpp
 * @brief Python module kessc: the SystemC model of key_equation_solver, cycle by cycle.
 *
 * SystemC elaborates once per process, so a single bench holds every configuration (two codes
 * times three architectures) driven by one clock, and is created on first use. Each call to
 * clock() is one rising edge with the same meaning as kes_timing.KesTimed.clock(): the given
 * inputs go to one instance and every other instance sees an idle edge.
 */

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "kes_sc.hpp"

namespace py = pybind11;
using namespace fec;

// libsystemc refers to sc_main; pybind11 builds with hidden visibility, so it must be exported.
__attribute__((visibility("default"))) int sc_main(int, char*[]) { return 0; }

namespace {

py::object bus_to_int(const std::string& bits) {
    return py::module_::import("builtins").attr("int")(bits, 2);
}

py::dict snapshot_to_dict(const Snapshot& s) {
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

/** @brief Type-erased access to one KesModule and its signals. */
class Instance {
public:
    virtual ~Instance() = default;
    virtual void drive(bool rst, bool valid, const Poly& syndromes) = 0;
    virtual py::dict ports() const = 0;
    virtual py::dict snapshot() const = 0;
    virtual py::dict config() const = 0;
};

template <unsigned M, unsigned T>
class Bench : public Instance {
public:
    using Dut = KesModule<M, T>;

    Bench(const std::string& name, const std::string& architecture, std::uint32_t poly,
          sc_core::sc_clock& clk)
        : field_(M, poly),
          architecture_(architecture),
          rst_((name + "_rst").c_str()),
          valid_((name + "_valid").c_str()),
          syn_((name + "_syn").c_str()),
          o_valid_((name + "_o_valid").c_str()),
          o_lambda_((name + "_o_lambda").c_str()),
          o_omega_((name + "_o_omega").c_str()),
          dut_(name.c_str(), architecture, field_) {
        dut_.i_clk(clk);
        dut_.i_rst(rst_);
        dut_.i_valid(valid_);
        dut_.i_syn(syn_);
        dut_.o_valid(o_valid_);
        dut_.o_lambda(o_lambda_);
        dut_.o_omega(o_omega_);
    }

    void drive(bool rst, bool valid, const Poly& syndromes) override {
        if (valid) {
            if (syndromes.size() != 2 * T) {
                throw std::invalid_argument("expected " + std::to_string(2 * T) + " syndromes, got " +
                                            std::to_string(syndromes.size()));
            }
            for (Symbol s : syndromes) {
                field_.check(s);
            }
            syn_.write(pack_bus<M, Dut::SYN_W>(syndromes));
        }
        rst_.write(rst);
        valid_.write(valid);
    }

    py::dict ports() const override {
        py::dict d;
        d["o_valid"] = o_valid_.read();
        d["o_lambda"] = py::cast(unpack_bus<M>(o_lambda_.read(), T + 1));
        d["o_omega"] = py::cast(unpack_bus<M>(o_omega_.read(), T));
        d["lambda_bus"] = bus_to_int(o_lambda_.read().to_string());
        d["omega_bus"] = bus_to_int(o_omega_.read().to_string());
        d["busy"] = dut_.busy();
        const auto state = dut_.state();
        d["state"] = state ? py::cast(*state) : py::none();
        return d;
    }

    py::dict snapshot() const override { return snapshot_to_dict(dut_.core().snapshot()); }

    py::dict config() const override {
        py::dict d;
        d["m"] = M;
        d["poly"] = field_.poly();
        d["t"] = T;
        d["architecture"] = architecture_;
        d["iterations"] = dut_.core().iterations();
        d["widths"] = py::make_tuple(Dut::SYN_W, Dut::LAMBDA_W, Dut::OMEGA_W);
        return d;
    }

private:
    GaloisField field_;
    std::string architecture_;
    sc_core::sc_signal<bool> rst_, valid_;
    sc_core::sc_signal<sc_dt::sc_bv<Dut::SYN_W>> syn_;
    sc_core::sc_signal<bool> o_valid_;
    sc_core::sc_signal<sc_dt::sc_bv<Dut::LAMBDA_W>> o_lambda_;
    sc_core::sc_signal<sc_dt::sc_bv<Dut::OMEGA_W>> o_omega_;
    Dut dut_;
};

/**
 * @brief Every configuration on one clock: 10 ns period, first rising edge at 5 ns, so inputs
 *        written between two calls are stable when the edge arrives.
 */
class Harness {
public:
    static Harness& get() {
        static Harness* harness = new Harness();    // never destroyed: SystemC dislikes teardown at exit
        return *harness;
    }

    py::dict clock(const std::string& code, const std::string& architecture, bool rst, bool valid,
                   const std::optional<Poly>& syndromes) {
        Instance& target = find(code, architecture);
        if (valid && !syndromes) {
            throw std::invalid_argument("i_valid is set but no syndromes were given");
        }
        target.drive(rst, valid, syndromes.value_or(Poly{}));    // validates before touching others
        for (auto& [key, inst] : instances_) {
            if (inst.get() != &target) {
                inst->drive(false, false, Poly{});
            }
        }
        sc_core::sc_start(period_);
        ++cycles_;
        return target.ports();
    }

    py::dict ports(const std::string& code, const std::string& architecture) {
        return find(code, architecture).ports();
    }
    py::dict snapshot(const std::string& code, const std::string& architecture) {
        return find(code, architecture).snapshot();
    }
    py::dict config(const std::string& code, const std::string& architecture) {
        return find(code, architecture).config();
    }

    std::vector<std::pair<std::string, std::string>> instances() const {
        std::vector<std::pair<std::string, std::string>> out;
        for (const auto& [key, inst] : instances_) {
            out.push_back(key);
        }
        return out;
    }

    double time_ns() const { return sc_core::sc_time_stamp() / sc_core::sc_time(1, sc_core::SC_NS); }
    unsigned long long cycles() const { return cycles_; }

private:
    Harness()
        : period_(10, sc_core::SC_NS),
          clk_("clk", period_, 0.5, sc_core::sc_time(5, sc_core::SC_NS), true) {
        for (const char* arch : {"ibm", "ribm", "ribm_systolic"}) {
            const std::string a(arch);
            instances_[{"qr", a}] = std::make_unique<Bench<8, 8>>("qr_" + a, a, 0x1D, clk_);
            instances_[{"tp3", a}] = std::make_unique<Bench<4, 3>>("tp3_" + a, a, 0x3, clk_);
        }
    }

    Instance& find(const std::string& code, const std::string& architecture) {
        auto it = instances_.find({code, architecture});
        if (it == instances_.end()) {
            throw std::invalid_argument("no instance '" + code + "', '" + architecture + "'");
        }
        return *it->second;
    }

    sc_core::sc_time period_;
    sc_core::sc_clock clk_;
    std::map<std::pair<std::string, std::string>, std::unique_ptr<Instance>> instances_;
    unsigned long long cycles_ = 0;
};

}  // namespace

PYBIND11_MODULE(kessc, m) {
    m.doc() = "SystemC model of key_equation_solver (see docs/kes_protocol.md)";

    py::class_<Harness, std::unique_ptr<Harness, py::nodelete>>(m, "Harness")
        .def("clock", &Harness::clock, py::arg("code"), py::arg("architecture"),
             py::arg("rst") = false, py::arg("valid") = false, py::arg("syndromes") = py::none(),
             "One rising edge: the inputs go to one instance, the others see an idle edge.")
        .def("ports", &Harness::ports, py::arg("code"), py::arg("architecture"))
        .def("snapshot", &Harness::snapshot, py::arg("code"), py::arg("architecture"))
        .def("config", &Harness::config, py::arg("code"), py::arg("architecture"))
        .def_property_readonly("instances", &Harness::instances)
        .def_property_readonly("time_ns", &Harness::time_ns)
        .def_property_readonly("cycles", &Harness::cycles);

    m.def("harness", &Harness::get, py::return_value_policy::reference,
          "The bench, elaborated on the first call.");
    m.def("release", [] { return std::string(sc_core::sc_release()); });
}
