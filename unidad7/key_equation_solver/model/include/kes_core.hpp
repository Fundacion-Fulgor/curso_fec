/**
 * @file kes_core.hpp
 * @brief Register-level models of the key equation solvers of Sarwate and Shanbhag (2001).
 *
 * Line by line translation of unidad7/kes_ref.py, with the same states, register names and
 * results, so both can be compared state by state:
 *
 *   IBM           inversionless BM, Blahut version, 3t steps, classic Omega
 *   RiBMSplit     riBM, rDC block + ELU block, 2t steps, high-order Omega
 *   RiBMSystolic  RiBM, one array of 3t + 1 PE1, 2t steps, high-order Omega
 *
 * load() sets state 0 from S_0 .. S_(2t-1), each step() is one clock edge, and result() reads the
 * outputs after the last step. Coefficient vectors are in ascending order.
 */
#pragma once

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gf.hpp"

namespace fec {

using Poly = std::vector<Symbol>;

/** @brief Register contents and control signals of one state. */
struct Snapshot {
    unsigned r = 0;
    Symbol gamma = 0;
    int k = 0;
    bool mc = false;
    Symbol discrepancy = 0;
    std::vector<std::pair<std::string, Poly>> registers;
};

/** @brief Output of a solver; see kes_ref.KesResult. */
struct KesResult {
    Poly lam;
    Poly omega;
    std::string form;          ///< "classic" or "high"
    unsigned iterations = 0;
    int length = 0;            ///< LFSR length L = (2t - k) / 2
    bool traced = false;
    std::vector<Snapshot> trace;

    /** @brief Degree of Lambda, -1 if it is zero. */
    int degree() const {
        for (int i = static_cast<int>(lam.size()) - 1; i >= 0; --i) {
            if (lam[static_cast<std::size_t>(i)]) {
                return i;
            }
        }
        return -1;
    }
};

/** @brief State handling shared by the three architectures. */
class Kes {
public:
    /**
     * @param field Symbol field (copied).
     * @param t Error correcting capability; the solver takes 2t syndromes.
     * @throws std::invalid_argument If t < 1.
     */
    Kes(const GaloisField& field, unsigned t) : field_(field), t_(t) {
        if (t < 1) {
            throw std::invalid_argument("t must be >= 1");
        }
    }
    virtual ~Kes() = default;

    virtual std::string name() const = 0;
    virtual std::string form() const = 0;
    virtual Symbol discrepancy() const = 0;
    /** @brief Clock edges from load() to the result. */
    virtual unsigned iterations() const { return 2 * t_; }

    const GaloisField& field() const { return field_; }
    unsigned t() const { return t_; }
    unsigned r() const { return r_; }
    bool loaded() const { return loaded_; }
    bool done() const { return loaded_ && r_ == iterations(); }
    Symbol gamma() const { return gamma_; }
    int k() const { return k_; }

    /** @brief MC(r): discrepancy != 0 and k >= 0, only during the 2t BM steps. */
    bool mc() const { return loaded_ && r_ < 2 * t_ && discrepancy() != 0 && k_ >= 0; }

    /** @brief Clear every register, as a synchronous reset does. */
    void reset() {
        r_ = 0;
        loaded_ = false;
        gamma_ = 0;
        k_ = 0;
        clear();
    }

    /**
     * @brief Load the syndromes (state 0).
     * @throws std::invalid_argument If there are not 2t syndromes or one is not a field element.
     */
    void load(const Poly& syndromes) {
        if (syndromes.size() != 2 * t_) {
            throw std::invalid_argument("expected " + std::to_string(2 * t_) + " syndromes, got " +
                                        std::to_string(syndromes.size()));
        }
        for (Symbol s : syndromes) {
            field_.check(s);
        }
        r_ = 0;
        loaded_ = true;
        gamma_ = 1;
        k_ = 0;
        init(syndromes);
    }

    /**
     * @brief One clock edge.
     * @throws std::runtime_error If nothing is loaded or every step is done.
     */
    void step() {
        if (!loaded_) {
            throw std::runtime_error("load() must be called before step()");
        }
        if (r_ >= iterations()) {
            throw std::runtime_error("all " + std::to_string(iterations()) + " steps are done");
        }
        advance();
        ++r_;
    }

    /**
     * @brief Outputs after the last step.
     * @throws std::runtime_error If the computation is not finished.
     */
    KesResult result() const {
        if (!done()) {
            throw std::runtime_error("the computation is not finished");
        }
        KesResult res;
        outputs(res.lam, res.omega);
        res.form = form();
        res.iterations = iterations();
        res.length = (2 * static_cast<int>(t_) - k_) / 2;
        return res;
    }

    /** @brief load() followed by every step, optionally keeping every state. */
    KesResult run(const Poly& syndromes, bool trace = false) {
        load(syndromes);
        std::vector<Snapshot> states;
        if (trace) {
            states.push_back(snapshot());
        }
        while (r_ < iterations()) {
            step();
            if (trace) {
                states.push_back(snapshot());
            }
        }
        KesResult res = result();
        res.traced = trace;
        res.trace = std::move(states);
        return res;
    }

    Snapshot snapshot() const {
        Snapshot s;
        s.r = r_;
        s.gamma = gamma_;
        s.k = k_;
        s.mc = mc();
        s.discrepancy = discrepancy();
        registers(s.registers);
        return s;
    }

protected:
    virtual void clear() = 0;
    virtual void init(const Poly& syndromes) = 0;
    virtual void advance() = 0;
    virtual void outputs(Poly& lam, Poly& omega) const = 0;
    virtual void registers(std::vector<std::pair<std::string, Poly>>& regs) const = 0;

    /** @brief Control unit: update gamma and k, return MC. */
    bool control(Symbol discrepancy) {
        const bool mc = discrepancy != 0 && k_ >= 0;
        if (mc) {
            gamma_ = discrepancy;
            k_ = -k_ - 1;
        } else {
            ++k_;
        }
        return mc;
    }

    /** @brief Datapath of PE0 and PE1: gamma * upper + d0 * lower. */
    Poly pe_update(const Poly& upper, const Poly& lower, Symbol d0) const {
        Poly out(upper.size());
        for (std::size_t i = 0; i < upper.size(); ++i) {
            out[i] = field_.mul(gamma_, upper[i]) ^ field_.mul(d0, lower[i]);
        }
        return out;
    }

    /** @brief p shifted one place towards the high end: [0, p_0, .., p_(n-2)]. */
    static Poly shift_up(const Poly& p) {
        Poly out(p.size(), 0);
        std::copy(p.begin(), p.end() - 1, out.begin() + 1);
        return out;
    }

    /** @brief p shifted one place towards the low end: [p_1, .., p_(n-1), 0]. */
    static Poly shift_down(const Poly& p) {
        Poly out(p.size(), 0);
        std::copy(p.begin() + 1, p.end(), out.begin());
        return out;
    }

    GaloisField field_;
    unsigned t_;
    unsigned r_ = 0;
    bool loaded_ = false;
    Symbol gamma_ = 0;
    int k_ = 0;
};

/**
 * @brief iBM, Blahut version (section III-B): DC block with a circular syndrome register DS and
 *        t + 1 multipliers, ELU block of t + 1 PE0. Steps 2t .. 3t-1 produce Omega.
 */
class IBM : public Kes {
public:
    IBM(const GaloisField& field, unsigned t) : Kes(field, t) { reset(); }

    std::string name() const override { return "ibm"; }
    std::string form() const override { return "classic"; }
    unsigned iterations() const override { return 3 * t_; }

    /** @brief DC block output: sum of DS_((2t - i) mod 2t) * lambda_i. */
    Symbol discrepancy() const override {
        const std::size_t n = 2 * t_;
        Symbol acc = 0;
        for (std::size_t i = 0; i < lam_.size(); ++i) {
            acc ^= field_.mul(ds_[(n - i % n) % n], lam_[i]);
        }
        return acc;
    }

protected:
    void clear() override {
        lam_.assign(t_ + 1, 0);
        b_.assign(t_ + 1, 0);
        ds_.assign(2 * t_, 0);
        omega_.clear();
    }

    void init(const Poly& syndromes) override {
        lam_.assign(t_ + 1, 0);
        lam_[0] = 1;
        b_ = lam_;
        ds_ = syndromes;
        omega_.clear();
    }

    void advance() override {
        const std::size_t n = 2 * t_;
        const Symbol delta = discrepancy();
        Poly rotated(ds_.begin() + 1, ds_.end());
        rotated.push_back(ds_[0]);
        if (r_ < n) {
            Poly lam = pe_update(lam_, shift_up(b_), delta);
            b_ = control(delta) ? lam_ : shift_up(b_);
            lam_ = std::move(lam);
            if (r_ == n - 1) {
                std::fill(rotated.begin() + t_, rotated.end(), 0);
            }
        } else {
            omega_.push_back(delta);
        }
        ds_ = std::move(rotated);
    }

    void outputs(Poly& lam, Poly& omega) const override {
        lam = lam_;
        omega = omega_;
    }

    void registers(std::vector<std::pair<std::string, Poly>>& regs) const override {
        regs = {{"lam", lam_}, {"b", b_}, {"ds", ds_}, {"omega", omega_}};
    }

private:
    Poly lam_, b_, ds_, omega_;
};

/** @brief riBM (section IV-B, Fig. 5): rDC block of 2t PE1 and ELU block of t + 1 PE0. */
class RiBMSplit : public Kes {
public:
    RiBMSplit(const GaloisField& field, unsigned t) : Kes(field, t) { reset(); }

    std::string name() const override { return "ribm"; }
    std::string form() const override { return "high"; }
    Symbol discrepancy() const override { return delta_[0]; }

protected:
    void clear() override {
        lam_.assign(t_ + 1, 0);
        b_.assign(t_ + 1, 0);
        delta_.assign(2 * t_, 0);
        theta_.assign(2 * t_, 0);
    }

    void init(const Poly& syndromes) override {
        lam_.assign(t_ + 1, 0);
        lam_[0] = 1;
        b_ = lam_;
        delta_ = syndromes;
        theta_ = syndromes;
    }

    void advance() override {
        const Symbol d0 = delta_[0];
        const Poly shifted = shift_down(delta_);
        Poly lam = pe_update(lam_, shift_up(b_), d0);
        Poly delta = pe_update(shifted, theta_, d0);
        if (control(d0)) {
            b_ = lam_;
            theta_ = shifted;
        } else {
            b_ = shift_up(b_);
        }
        lam_ = std::move(lam);
        delta_ = std::move(delta);
    }

    void outputs(Poly& lam, Poly& omega) const override {
        lam = lam_;
        omega.assign(delta_.begin(), delta_.begin() + t_);
    }

    void registers(std::vector<std::pair<std::string, Poly>>& regs) const override {
        regs = {{"lam", lam_}, {"b", b_}, {"delta", delta_}, {"theta", theta_}};
    }

private:
    Poly lam_, b_, delta_, theta_;
};

/** @brief RiBM (section IV-B3, Fig. 6): a single array of 3t + 1 PE1. */
class RiBMSystolic : public Kes {
public:
    RiBMSystolic(const GaloisField& field, unsigned t) : Kes(field, t) { reset(); }

    std::string name() const override { return "ribm_systolic"; }
    std::string form() const override { return "high"; }
    Symbol discrepancy() const override { return delta_[0]; }

protected:
    void clear() override {
        delta_.assign(3 * t_ + 1, 0);
        theta_.assign(3 * t_ + 1, 0);
    }

    void init(const Poly& syndromes) override {
        delta_ = syndromes;
        delta_.resize(3 * t_ + 1, 0);
        delta_.back() = 1;
        theta_ = delta_;
    }

    void advance() override {
        const Symbol d0 = delta_[0];
        const Poly shifted = shift_down(delta_);
        delta_ = pe_update(shifted, theta_, d0);
        if (control(d0)) {
            theta_ = shifted;
        }
    }

    void outputs(Poly& lam, Poly& omega) const override {
        lam.assign(delta_.begin() + t_, delta_.begin() + 2 * t_ + 1);
        omega.assign(delta_.begin(), delta_.begin() + t_);
    }

    void registers(std::vector<std::pair<std::string, Poly>>& regs) const override {
        regs = {{"delta", delta_}, {"theta", theta_}};
    }

private:
    Poly delta_, theta_;
};

/**
 * @brief Build a solver by name.
 * @param name "ibm", "ribm" or "ribm_systolic".
 * @throws std::invalid_argument If the name is unknown.
 */
inline std::unique_ptr<Kes> make_kes(const std::string& name, const GaloisField& field, unsigned t) {
    if (name == "ibm") {
        return std::make_unique<IBM>(field, t);
    }
    if (name == "ribm") {
        return std::make_unique<RiBMSplit>(field, t);
    }
    if (name == "ribm_systolic") {
        return std::make_unique<RiBMSystolic>(field, t);
    }
    throw std::invalid_argument("unknown architecture '" + name +
                                "', expected one of ['ibm', 'ribm', 'ribm_systolic']");
}

}  // namespace fec
