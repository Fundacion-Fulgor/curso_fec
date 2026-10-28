/**
 * @file kes_sc.hpp
 * @brief SystemC model of key_equation_solver: the ports of the RTL and the protocol of
 *        docs/kes_protocol.md around the register-level solvers of kes_core.hpp.
 *
 * Buses hold the degree-i coefficient in bits [M*i + M - 1 : M*i]. With M = 8 and T = 8 the
 * widths are those of TP3: i_syn 128, o_lambda 72, o_omega 64.
 */

#pragma once

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include <systemc>

#include "kes_core.hpp"

namespace fec {

/** @brief Unpack count M-bit symbols from a bus. */
template <unsigned M, int W>
Poly unpack_bus(const sc_dt::sc_bv<W>& bus, unsigned count) {
    Poly out(count);
    for (unsigned i = 0; i < count; ++i) {
        out[i] = static_cast<Symbol>(bus.range(M * i + M - 1, M * i).to_uint());
    }
    return out;
}

/** @brief Pack M-bit symbols into a bus, symbol i in bits [M*i + M - 1 : M*i]. */
template <unsigned M, int W>
sc_dt::sc_bv<W> pack_bus(const Poly& symbols) {
    sc_dt::sc_bv<W> bus(0);
    for (unsigned i = 0; i < symbols.size(); ++i) {
        bus.range(M * i + M - 1, M * i) = symbols[i];
    }
    return bus;
}

/**
 * @brief key_equation_solver for GF(2^M) and capability T.
 *
 * One SC_METHOD on the rising edge of i_clk applies the rules of the protocol in order:
 * reset, one step of the word in progress (with the result and o_valid after the last step),
 * and the load of a new word when i_valid is high.
 */
template <unsigned M, unsigned T>
class KesModule : public sc_core::sc_module {
public:
    static constexpr int SYN_W = static_cast<int>(2 * T * M);
    static constexpr int LAMBDA_W = static_cast<int>((T + 1) * M);
    static constexpr int OMEGA_W = static_cast<int>(T * M);

    sc_core::sc_in<bool> i_clk{"i_clk"};
    sc_core::sc_in<bool> i_rst{"i_rst"};
    sc_core::sc_in<bool> i_valid{"i_valid"};
    sc_core::sc_in<sc_dt::sc_bv<SYN_W>> i_syn{"i_syn"};
    sc_core::sc_out<bool> o_valid{"o_valid"};
    sc_core::sc_out<sc_dt::sc_bv<LAMBDA_W>> o_lambda{"o_lambda"};
    sc_core::sc_out<sc_dt::sc_bv<OMEGA_W>> o_omega{"o_omega"};

    /**
     * @param name Module name.
     * @param architecture "ibm", "ribm" or "ribm_systolic".
     * @param field Symbol field; its degree must be M.
     * @throws std::invalid_argument If the field degree is not M or the architecture is unknown.
     */
    KesModule(sc_core::sc_module_name name, const std::string& architecture, const GaloisField& field)
        : sc_core::sc_module(name), core_(make_kes(architecture, field, T)) {
        if (field.m() != M) {
            throw std::invalid_argument("the field degree does not match the module");
        }
        SC_METHOD(on_clock);
        sensitive << i_clk.pos();
        dont_initialize();
    }

    const Kes& core() const { return *core_; }
    bool busy() const { return busy_; }

    /** @brief Index of the internal state while a word is in progress. */
    std::optional<unsigned> state() const {
        return busy_ ? std::optional<unsigned>(core_->r()) : std::nullopt;
    }

private:
    void on_clock() {
        if (i_rst.read()) {
            core_->reset();
            busy_ = false;
            o_valid.write(false);
            o_lambda.write(sc_dt::sc_bv<LAMBDA_W>(0));
            o_omega.write(sc_dt::sc_bv<OMEGA_W>(0));
            return;
        }
        bool valid = false;
        if (busy_) {
            core_->step();
            if (core_->done()) {
                const KesResult res = core_->result();
                o_lambda.write(pack_bus<M, LAMBDA_W>(res.lam));
                o_omega.write(pack_bus<M, OMEGA_W>(res.omega));
                valid = true;
                busy_ = false;
            }
        }
        if (i_valid.read()) {
            core_->load(unpack_bus<M>(i_syn.read(), 2 * T));
            busy_ = true;
        }
        o_valid.write(valid);
    }

    std::unique_ptr<Kes> core_;
    bool busy_ = false;
};

}  // namespace fec
