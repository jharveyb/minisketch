/**********************************************************************
 * Copyright (c) 2018 Pieter Wuille, Greg Maxwell, Gleb Naumenko      *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include "../include/minisketch.h"
#include <string.h>
#include <memory>
#include <vector>
#include <chrono>
#include <random>
#include <set>
#include <algorithm>

// Human-readable name for each implementation ID; must stay in sync with
// FieldImpl in src/minisketch.cpp.
static const char* ImplName(uint32_t impl) {
    switch (impl) {
        case 0: return "GENERIC";
        case 1: return "CLMUL";
        case 2: return "CLMUL_TRI";
        default: return "UNKNOWN";
    }
}

// Print CLI usage; invoked on --help/-h or on any argument error.
static void PrintUsage(const char* prog) {
    printf("Usage: %s [--sweep] [--sweep-errors] [--create-only] [--data-sweep] [--impl N] [--data N] [--syndromes N] [--errors N] [--iters N]\n", prog);
    printf("  --sweep          run fixed matrix: bits in {64}, capacity in {1024,2048,...,8192}\n");
    printf("                   (in sweep mode --syndromes is ignored)\n");
    printf("  --sweep-errors   for each swept capacity, test errors at 10%%, 25%%, 50%%, 75%%, 100%% of capacity\n");
    printf("                   (implies --sweep; ignores --errors)\n");
    printf("  --create-only    skip the 'recover' benchmark, only run the 'create' benchmark\n");
    printf("  --data-sweep     sweep data_len over {8192, 16384, 51200, 102400, 256000, 512000}\n");
    printf("                   (implies --create-only; conflicts with --data and --sweep-errors)\n");
    printf("  --impl N         restrict to implementation N (0..minisketch_implementation_max())\n");
    printf("                   default: iterate all implementations\n");
    printf("  --data N         elements added per sketch in the 'create' phase (default: errors*10)\n");
    printf("  --syndromes N    sketch capacity in non-sweep mode (default: 150, range 0..1000000)\n");
    printf("  --errors N       number of elements per sketch (default: syndromes; range >= 0)\n");
    printf("  --iters N        benchmark iterations (default: 10, range 0..1000000000)\n");
}

// Emit the column header exactly once so stdout is a valid long-format TSV.
static void EmitTsvHeader() {
    printf("metric\tbits\tcapacity\terrors\tdata_len\timplementation\tvalue\n");
}

// Benchmark minisketch_decode for each impl in [impl_lo, impl_hi]. Emits one
// recover[ms] row per impl with the best (minimum) time across `iters` runs.
// Impls unsupported for this (bits, CPU) are silently skipped.
static void RunDecodeBench(int bits, int syndromes, int errors, int iters,
                           uint32_t impl_lo, uint32_t impl_hi) {
    for (uint32_t impl = impl_lo; impl <= impl_hi; ++impl) {
        std::vector<minisketch*> states;
        std::vector<uint64_t> roots(2 * syndromes);
        std::random_device rng;
        std::uniform_int_distribution<uint64_t> dist(1, (uint64_t(1) << bits) - 1);
        states.resize(iters);
        std::vector<double> benches;
        benches.reserve(iters);
        for (int i = 0; i < iters; ++i) {
            states[i] = minisketch_create(bits, impl, syndromes);
            if (!states[i]) break;
            std::set<uint64_t> done;
            for (int j = 0; j < errors; ++j) {
                uint64_t r;
                do {
                    r = dist(rng);
                } while (done.count(r));
                done.insert(r);
                minisketch_add_uint64(states[i], r);
            }
        }
        if (states[0]) {
            for (auto& state : states) {
                auto start = std::chrono::steady_clock::now();
                minisketch_decode(state, 2 * syndromes, roots.data());
                auto stop = std::chrono::steady_clock::now();
                std::chrono::duration<double> dur(stop - start);
                benches.push_back(dur.count());
            }
            std::sort(benches.begin(), benches.end());
            printf("recover[ms]\t%i\t%i\t%i\t0\t%s\t%.5f\n",
                   bits, syndromes, errors, ImplName(impl), benches[0] * 1000.0);
        }
        for (auto& state : states) {
            minisketch_destroy(state);
        }
    }
}

// Benchmark minisketch_add_uint64 for each impl in [impl_lo, impl_hi]. Emits
// two rows per impl: create[ns] normalises best time by (data_len * syndromes)
// to expose the per-syndrome primitive cost, and create-total[ms] reports the
// absolute wall time for all add calls.
static void RunCreateBench(int bits, int syndromes, int errors, int iters,
                           long data_len, uint32_t impl_lo, uint32_t impl_hi) {
    for (uint32_t impl = impl_lo; impl <= impl_hi; ++impl) {
        std::vector<minisketch*> states;
        std::random_device rng;
        std::uniform_int_distribution<uint64_t> dist;
        std::vector<uint64_t> data;
        data.resize(data_len);
        states.resize(iters);
        std::vector<double> benches;
        benches.reserve(iters);
        for (int i = 0; i < iters; ++i) {
            states[i] = minisketch_create(bits, impl, syndromes);
        }
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] = dist(rng);
        }
        if (states[0]) {
            for (auto& state : states) {
                auto start = std::chrono::steady_clock::now();
                for (auto val : data) {
                    minisketch_add_uint64(state, val);
                }
                auto stop = std::chrono::steady_clock::now();
                std::chrono::duration<double> dur(stop - start);
                benches.push_back(dur.count());
            }
            std::sort(benches.begin(), benches.end());
            printf("create[ns]\t%i\t%i\t%i\t%ld\t%s\t%.5f\n",
                   bits, syndromes, errors, data_len, ImplName(impl),
                   benches[0] * 1000000000.0 / data_len / syndromes);
            printf("create-total[ms]\t%i\t%i\t%i\t%ld\t%s\t%.5f\n",
                   bits, syndromes, errors, data_len, ImplName(impl),
                   benches[0] * 1000.0);
        }
        for (auto& state : states) {
            minisketch_destroy(state);
        }
    }
}

// Accept either "--flag value" (consuming argv[++i]) or "--flag=value".
// Returns true and sets out_val on match.
static bool ParseFlagValue(const char* arg, const char* prefix, size_t prefix_len,
                           int argc, char** argv, int& i, const char*& out_val) {
    if (strcmp(arg, prefix) == 0) {
        if (i + 1 >= argc) return false;
        out_val = argv[++i];
        return true;
    }
    if (strncmp(arg, prefix, prefix_len) == 0 && arg[prefix_len] == '=') {
        out_val = arg + prefix_len + 1;
        return true;
    }
    return false;
}

// Parse flags, validate config, emit the TSV header, then dispatch to the
// requested sweep/non-sweep and decode/create combinations.
int main(int argc, char** argv) {
    bool sweep = false;
    bool sweep_errors = false;
    bool create_only = false;
    bool data_sweep = false;
    bool explicit_bits = false;
    int impl_pin = -1;
    long data_override = -1;
    int syndromes = 150;
    int errors = -1;
    int iters = 10;
    int bits = 64;

    for (int i = 1; i < argc; ++i) {
        char* a = argv[i];
        const char* v = nullptr;
        if (strcmp(a, "--help") == 0 || strcmp(a, "-h") == 0) {
            PrintUsage(argv[0]);
            return 0;
        } else if (strcmp(a, "--sweep") == 0) {
            sweep = true;
        } else if (strcmp(a, "--sweep-errors") == 0) {
            sweep_errors = true;
            sweep = true;
        } else if (strcmp(a, "--create-only") == 0) {
            create_only = true;
        } else if (strcmp(a, "--data-sweep") == 0) {
            data_sweep = true;
            create_only = true;
        } else if (ParseFlagValue(a, "--impl", 6, argc, argv, i, v)) {
            impl_pin = strtol(v, NULL, 10);
        } else if (ParseFlagValue(a, "--data", 6, argc, argv, i, v)) {
            data_override = strtol(v, NULL, 10);
        } else if (ParseFlagValue(a, "--syndromes", 11, argc, argv, i, v)) {
            syndromes = (int)strtol(v, NULL, 10);
        } else if (ParseFlagValue(a, "--errors", 8, argc, argv, i, v)) {
            errors = (int)strtol(v, NULL, 10);
        } else if (ParseFlagValue(a, "--iters", 7, argc, argv, i, v)) {
            iters = (int)strtol(v, NULL, 10);
        } else if (ParseFlagValue(a, "--bits", 6, argc, argv, i, v)) {
            explicit_bits = true;
            bits = (int)strtol(v, NULL, 10);
        } else {
            printf("Unexpected argument: %s\n", a);
            PrintUsage(argv[0]);
            return 1;
        }
    }

    bool errors_given = (errors >= 0);
    if (!errors_given) errors = syndromes;

    if (!sweep && (syndromes < 0 || syndromes > 1000000)) {
        printf("--syndromes (%i) out of range 0..1000000\n", syndromes);
        return 1;
    }
    if (errors_given && errors < 0) {
        printf("--errors (%i) is negative\n", errors);
        return 1;
    }
    if (iters < 0 || iters > 1000000000) {
        printf("--iters (%i) out of range 0..1000000000\n", iters);
        return 1;
    }
    if (data_override != -1 && (data_override < 1 || data_override > 1000000000)) {
        printf("--data value (%ld) out of range 1..1000000000\n", data_override);
        return 1;
    }
    if (data_sweep && data_override != -1) {
        printf("--data-sweep conflicts with --data\n");
        return 1;
    }
    if (data_sweep && sweep_errors) {
        printf("--data-sweep conflicts with --sweep-errors\n");
        return 1;
    }

    uint32_t max_impl = minisketch_implementation_max();
    if (impl_pin != -1) {
        if (impl_pin < 0 || (uint32_t)impl_pin > max_impl) {
            printf("--impl value (%d) out of range 0..%u\n", impl_pin, max_impl);
            return 1;
        }
    }

    uint32_t impl_lo = (impl_pin >= 0) ? (uint32_t)impl_pin : 0u;
    uint32_t impl_hi = (impl_pin >= 0) ? (uint32_t)impl_pin : max_impl;

    const long data_lens_sweep[] = {16, 64, 128, 512, 1024, 2048, 4096, 8192, 16384, 51200, 102400};

    if (sweep_errors && errors_given) {
        fprintf(stderr, "Note: in --sweep-errors mode, --errors is ignored\n");
    }
    if (data_sweep && errors_given) {
        fprintf(stderr, "Note: in --data-sweep mode, --errors is ignored\n");
    }

    EmitTsvHeader();

    // Covers multiple branches of the CLMUL implementation.
    // 2 bytes, 4 bytes with CLMUL_TRI support, 4 bytes with CLMUL where BITS <= 32,
    // 6 bytes with CLMUL where BITS % 8 == 0, 8 bytes with CLMUL_TRI, and finally
    // BITS == 64 for CLMUL.
    std::vector<int> default_bits_sweep = {16, 29, 32, 48, 58, 64};
    std::vector<int> custom_bits_sweep = {bits, 0, 0, 0, 0, 0};
    if (sweep) {
        const int caps[] = {128, 256, 512, 1024, 2048, 3072, 4096, 5120, 6144, 7168, 8192};
        std::vector<int> const bits_sweep = explicit_bits ? custom_bits_sweep : default_bits_sweep;
        const int errcounts[] = {8, 16, 32, 64, 128, 256, 384, 512, 768, 1024, 1536, 2048};
        for (int bits : bits_sweep) {
            if (!minisketch_bits_supported(bits)) continue;
            for (int cap : caps) {
                if (data_sweep) {
                    for (long data_len : data_lens_sweep) {
                        RunCreateBench(bits, cap, /*errors=*/0, iters, data_len,
                                       impl_lo, impl_hi);
                    }
                    continue;
                }
                std::vector<int> err_list;
                if (sweep_errors) {
                    for (int errcount : errcounts) {
                        if (errcount <= cap) {
                            err_list.push_back(errcount);
                        } 
                    }
                } else {
                    err_list.push_back(errors_given ? errors : cap);
                }
                for (int eff_errors : err_list) {
                    if (eff_errors > cap) continue;
                    if (eff_errors > pow(2.0, bits - 1)) continue;
                    long data_len = (data_override >= 0) ? data_override : (long)eff_errors * 10;
                    if (!create_only) {
                        RunDecodeBench(bits, cap, eff_errors, iters, impl_lo, impl_hi);
                    }
                    RunCreateBench(bits, cap, eff_errors, iters, data_len, impl_lo, impl_hi);
                }
            }
        }
    } else {
        for (int bits = 2; bits <= 64; ++bits) {
            if (!minisketch_bits_supported(bits)) continue;
            if (data_sweep) {
                for (long data_len : data_lens_sweep) {
                    RunCreateBench(bits, syndromes, /*errors=*/0, iters, data_len,
                                   impl_lo, impl_hi);
                }
                continue;
            }
            if (errors > pow(2.0, bits - 1)) continue;
            long data_len = (data_override >= 0) ? data_override : (long)errors * 10;
            if (!create_only) {
                RunDecodeBench(bits, syndromes, errors, iters, impl_lo, impl_hi);
            }
            RunCreateBench(bits, syndromes, errors, iters, data_len, impl_lo, impl_hi);
        }
    }
    return 0;
}
