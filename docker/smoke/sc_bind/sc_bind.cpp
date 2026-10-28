// Smoke test: a clocked SystemC model driven from Python through pybind11.
#include <pybind11/pybind11.h>
#include <systemc>

// libsystemc refers to sc_main, so the module must define and export it (pybind11 builds
// with hidden visibility) even though Python never calls it. Without it, import fails with
// "undefined symbol: sc_main".
__attribute__((visibility("default"))) int sc_main(int, char*[]) { return 0; }

namespace {

SC_MODULE(Counter) {
    sc_core::sc_in<bool> clk;
    unsigned count = 0;

    SC_CTOR(Counter) {
        SC_METHOD(tick);
        sensitive << clk.pos();
        dont_initialize();
    }

    void tick() { ++count; }
};

// SystemC elaborates once per process, so the model is built on first use and kept.
class Harness {
public:
    static Harness& get() {
        static Harness h;
        return h;
    }

    unsigned run(unsigned cycles) {
        sc_core::sc_start(cycles * period_);
        return counter_.count;
    }

private:
    Harness() : period_(10, sc_core::SC_NS), clk_("clk", period_), counter_("counter") {
        counter_.clk(clk_);
    }

    sc_core::sc_time period_;
    sc_core::sc_clock clk_;
    Counter counter_;
};

}  // namespace

PYBIND11_MODULE(sc_bind, m) {
    m.doc() = "SystemC smoke test";
    m.def("run", [](unsigned cycles) { return Harness::get().run(cycles); },
          "Advance the clock and return the number of rising edges seen.");
    m.def("now_ns", [] { return sc_core::sc_time_stamp().to_seconds() * 1e9; });
    m.def("release", [] { return std::string(sc_core::sc_release()); });
}
