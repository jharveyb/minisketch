/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** oracle_dump: print the decoder's internal stage outputs for cases read
 * from stdin, for tests/conformance.py to compare against the Python oracle
 * in tests/pyminisketch.py.
 *
 * One case per input line, whitespace-separated decimal integers:
 *
 *     <op> <bits> <impl> <param> <na> a_1 .. a_na [<nb> b_1 .. b_nb]
 *
 * impl is 0 (generic), 1 (clmul) or 2 (clmul trinomial); all field elements
 * are canonical (as serialized / minisketch_add_uint64). Ops:
 *
 *   cantor  -                    the additive FFT's Cantor basis
 *   decode  a = odd syndromes,   BerlekampMasseyOdd, BerlekampMasseyRef and
 *           param = max_count,   the full SketchImpl decode
 *           b = [basis seed]
 *   roots   a = monic poly (low  orbit full-split flag, betas and traces;
 *           to high), param =    FindRoots without / with power sums
 *           basis, b = odd power
 *           sums (optional)
 *
 * Output: one JSON object per line (sorted root lists, null on failure,
 * "stats_*" = MINISKETCH_STAT counters of that stage, showing which tiers
 * ran). */

#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "fields/generic_1byte.cpp"
#include "fields/generic_2bytes.cpp"
#include "fields/generic_4bytes.cpp"

#include "oracle_dump_impl.h"

#if defined(HAVE_CLMUL) && (defined(__x86_64__) || defined(__i386__))
#include <cpuid.h>
#endif

namespace {
std::map<std::string, unsigned long> g_stats;
}

void MinisketchStageStat(const char* name) { ++g_stats[name]; }

namespace oracle_dump {

std::string TakeStats() {
    std::string s = "{";
    for (const auto& kv : g_stats) {
        if (s.size() > 1) s += ",";
        s += "\"" + kv.first + "\":" + std::to_string(kv.second);
    }
    g_stats.clear();
    return s + "}";
}

/** Defined in oracle_dump_clmul.cpp; "" if (bits, impl) is not available. */
std::string RunClmul(int impl, int bits, const Case& c);

std::string RunGeneric(int bits, const Case& c) {
    switch (bits) {
#ifdef ENABLE_FIELD_INT_2
    case 2: return Run<Field2>(0, 2, c);
#endif
#ifdef ENABLE_FIELD_INT_3
    case 3: return Run<Field3>(0, 3, c);
#endif
#ifdef ENABLE_FIELD_INT_4
    case 4: return Run<Field4>(0, 4, c);
#endif
#ifdef ENABLE_FIELD_INT_5
    case 5: return Run<Field5>(0, 5, c);
#endif
#ifdef ENABLE_FIELD_INT_6
    case 6: return Run<Field6>(0, 6, c);
#endif
#ifdef ENABLE_FIELD_INT_7
    case 7: return Run<Field7>(0, 7, c);
#endif
#ifdef ENABLE_FIELD_INT_8
    case 8: return Run<Field8>(0, 8, c);
#endif
#ifdef ENABLE_FIELD_INT_9
    case 9: return Run<Field9>(0, 9, c);
#endif
#ifdef ENABLE_FIELD_INT_10
    case 10: return Run<Field10>(0, 10, c);
#endif
#ifdef ENABLE_FIELD_INT_11
    case 11: return Run<Field11>(0, 11, c);
#endif
#ifdef ENABLE_FIELD_INT_12
    case 12: return Run<Field12>(0, 12, c);
#endif
#ifdef ENABLE_FIELD_INT_13
    case 13: return Run<Field13>(0, 13, c);
#endif
#ifdef ENABLE_FIELD_INT_14
    case 14: return Run<Field14>(0, 14, c);
#endif
#ifdef ENABLE_FIELD_INT_15
    case 15: return Run<Field15>(0, 15, c);
#endif
#ifdef ENABLE_FIELD_INT_16
    case 16: return Run<Field16>(0, 16, c);
#endif
#ifdef ENABLE_FIELD_INT_32
    case 32: return Run<Field32>(0, 32, c);
#endif
    }
    return "";
}

bool ClmulUsable() {
#if !defined(HAVE_CLMUL)
    return false;
#elif defined(__x86_64__) || defined(__i386__)
    unsigned eax, ebx, ecx, edx;
    return __get_cpuid(1, &eax, &ebx, &ecx, &edx) && (ecx & 0x2);
#else
    return true; // Trust the compile-time gate elsewhere.
#endif
}

} // namespace oracle_dump

int main() {
    using namespace oracle_dump;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        Case c;
        int bits = 0, impl = 0;
        size_t na = 0, nb = 0;
        if (!(in >> c.op)) continue;
        in >> bits >> impl >> c.param >> na;
        c.a.resize(na);
        for (auto& v : c.a) in >> v;
        if (in >> nb) {
            c.b.resize(nb);
            for (auto& v : c.b) in >> v;
        }
        if (!in.eof() && in.fail()) {
            std::cout << "{\"error\":\"parse\"}" << std::endl;
            continue;
        }
        std::string out;
        if (c.op == "roots" && (c.a.empty() || c.a.back() != 1 || c.param == 0)) {
            out = "{\"error\":\"roots needs a monic polynomial and a nonzero basis\"}";
        } else if (impl == 0) {
            out = RunGeneric(bits, c);
        } else if (ClmulUsable()) {
            out = RunClmul(impl, bits, c);
        }
        if (out.empty()) out = "{\"error\":\"unsupported\"}";
        std::cout << out << std::endl;
    }
    return 0;
}
