/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** The per-field part of oracle_dump (see oracle_dump.cpp): runs one case
 * through the decoder's internal stages for a field class F and prints the
 * stage outputs as one JSON line. Included from a TU that has included the
 * field implementation files (generic and clmul ones go in separate TUs, as
 * both define a `Field` template in the anonymous namespace). */

#ifndef _MINISKETCH_ORACLE_DUMP_IMPL_H_
#define _MINISKETCH_ORACLE_DUMP_IMPL_H_

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "sketch_impl.h"
#include "test_refimpl.h"

namespace oracle_dump {

/** One input line, already split: op, the op's scalar parameter, and the
 *  canonical field elements (see oracle_dump.cpp for the line format). */
struct Case {
    std::string op;
    uint64_t param = 0;
    std::vector<uint64_t> a, b;
};

/** Stage counters (MINISKETCH_STAT) since the last call, as a JSON object. */
std::string TakeStats();

inline std::string Json(const std::vector<uint64_t>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) s += ",";
        s += std::to_string(v[i]);
    }
    return s + "]";
}

/** A connection polynomial or root list; {} means failure (null), unless
 *  empty_ok (an empty root list is a valid result). */
inline std::string JsonOrNull(const std::vector<uint64_t>& v, bool empty_ok = false) {
    return (v.empty() && !empty_ok) ? "null" : Json(v);
}

template<typename F>
struct Runner {
    typedef typename F::Elem Elem;
    const F field;
    const int impl, bits;

    Runner(int impl_in, int bits_in) : field(), impl(impl_in), bits(bits_in) {}

    std::vector<Elem> In(const std::vector<uint64_t>& v) const {
        std::vector<Elem> out;
        for (uint64_t x : v) out.push_back(field.FromUint64(x));
        return out;
    }
    std::vector<uint64_t> Out(const std::vector<Elem>& v, bool sort = false) const {
        std::vector<uint64_t> out;
        for (Elem x : v) out.push_back(field.ToUint64(x));
        if (sort) std::sort(out.begin(), out.end());
        return out;
    }

    /** cantor: the additive FFT's Cantor basis (fields of 2^k bits). */
    std::string Cantor() const {
        if (bits & (bits - 1)) return "{\"basis\":null}";
        AdditiveFFT<F> fft;
        fft.Extend(field, bits);
        return "{\"basis\":" + Json(Out(fft.Basis())) + "}";
    }

    /** decode: odd syndromes a, max_degree param. BerlekampMasseyOdd (as decode
     *  calls it), the textbook BerlekampMasseyRef on all syndromes, and the
     *  full decode through SketchImpl (max_count param, basis seed b[0]). */
    std::string Decode(const Case& c) const {
        auto osyn = In(c.a);
        const size_t max_degree = c.param;
        TakeStats();
        auto conn = BerlekampMasseyOdd(osyn, 2 * osyn.size(), max_degree, field);
        std::string bm_stats = TakeStats();
        auto ref = ut::BerlekampMasseyRef(ReconstructAllSyndromes(osyn, field), max_degree, field);

        // The public decode, from the canonical serialization (as in
        // pyminisketch: element i at bit offset bits * i, little endian).
        SketchImpl<F> sketch(impl, bits);
        sketch.Init(osyn.size());
        std::vector<unsigned char> ser((osyn.size() * bits + 7) / 8 + 8, 0);
        for (size_t i = 0; i < c.a.size(); ++i) {
            for (int j = 0; j < bits; ++j) {
                if ((c.a[i] >> j) & 1) {
                    size_t pos = i * bits + j;
                    ser[pos / 8] |= (unsigned char)(1 << (pos % 8));
                }
            }
        }
        sketch.Deserialize(ser.data());
        sketch.SetSeed(c.b.empty() ? uint64_t(-1) : c.b[0]);
        std::vector<uint64_t> out(max_degree + 1);
        int n = sketch.Decode(max_degree, out.data());
        std::string decode_stats = TakeStats();
        std::string dec = "null";
        if (n >= 0) {
            out.resize(n);
            std::sort(out.begin(), out.end());
            dec = Json(out);
        }
        return "{\"bm\":" + JsonOrNull(Out(conn)) + ",\"bm_ref\":" + JsonOrNull(Out(ref)) +
               ",\"decode\":" + dec + ",\"stats_bm\":" + bm_stats + ",\"stats_decode\":" + decode_stats + "}";
    }

    /** roots: monic polynomial a (low to high), basis param, optional odd
     *  power sums b. The orbit's full-split flag with its betas and traces
     *  (ComputeTraceOrbit, degree >= 3), and the sorted FindRoots result
     *  without and (if b is given) with the power sums. */
    std::string Roots(const Case& c) const {
        auto poly = In(c.a);
        const Elem basis = field.FromUint64(c.param);
        std::string res = "{";
        const size_t d = poly.size() - 1;
        if (d >= 3) {
            size_t depths = std::min<size_t>(bits, 6);
            std::vector<std::vector<Elem>> traces;
            std::vector<Elem> betas;
            TakeStats();
            bool split = ComputeTraceOrbit(poly, basis, depths, field, traces, betas);
            res += std::string("\"split\":") + (split ? "true" : "false") + ",\"stats_orbit\":" + TakeStats();
            if (split) {
                res += ",\"betas\":" + Json(Out(betas)) + ",\"traces\":[";
                for (size_t k = 0; k < traces.size(); ++k) res += (k ? "," : "") + Json(Out(traces[k]));
                res += "]";
            }
            res += ",";
        }
        TakeStats();
        auto roots = FindRoots(poly, basis, field);
        res += "\"roots\":" + JsonOrNull(Out(roots, true), d == 0) + ",\"stats_roots\":" + TakeStats();
        if (!c.b.empty()) {
            auto psum = In(c.b);
            auto roots_ps = FindRoots(poly, basis, field, &psum);
            res += ",\"roots_psum\":" + JsonOrNull(Out(roots_ps, true), d == 0) + ",\"stats_roots_psum\":" + TakeStats();
        }
        return res + "}";
    }

    std::string Run(const Case& c) const {
        if (c.op == "cantor") return Cantor();
        if (c.op == "decode") return Decode(c);
        if (c.op == "roots") return Roots(c);
        return "{\"error\":\"unknown op\"}";
    }
};

template<typename F>
std::string Run(int impl, int bits, const Case& c) {
    Runner<F> r(impl, bits);
    return r.Run(c);
}

} // namespace oracle_dump

#endif
