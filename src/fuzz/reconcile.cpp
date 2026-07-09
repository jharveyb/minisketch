/**********************************************************************
 * Copyright (c) 2021-present The Bitcoin Core developers              *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Fuzz target: set reconciliation between two parties.
 *
 * Adapted from Bitcoin Core's src/test/fuzz/minisketch.cpp. Two sketches are
 * filled with a data-driven mix of shared and unique elements while the true
 * symmetric difference is tracked in a map; the merged difference sketch
 * (optionally after a serialization roundtrip on either side) must decode to
 * exactly that difference whenever it fits the capacity.
 *
 * Compared to the decode/roundtrip targets this uses only the two
 * practically-relevant field sizes (32 bits, like Bitcoin's
 * txreconciliation, and 64 bits) but larger capacities (log-uniform up to
 * 1024), exercising high-degree Berlekamp-Massey and root finding.
 */

#include <cstdint>
#include <map>
#include <numeric>
#include <utility>

#include "../../include/minisketch.h"
#include "FuzzedDataProvider.h"
#include "fuzz.h"

FUZZ_TARGET(reconcile) {
    FuzzedDataProvider fuzzed_data_provider{data, size};

    const uint32_t bits{fuzzed_data_provider.ConsumeBool() ? 64u : 32u};
    const uint64_t max_entry{bits == 64 ? std::numeric_limits<uint64_t>::max() - 1
                                        : std::numeric_limits<uint32_t>::max() - 1};
    // Log-uniform capacity in [1, 1024] (see decode.cpp).
    const size_t cap_lo{size_t{1} << fuzzed_data_provider.ConsumeIntegralInRange<int>(0, 10)};
    const auto capacity{fuzzed_data_provider.ConsumeIntegralInRange<size_t>(cap_lo, std::min<size_t>(2 * cap_lo - 1, 1024))};
    const uint32_t impl{fuzzed_data_provider.ConsumeIntegralInRange<uint32_t>(0, Minisketch::MaxImplementation())};
    if (!Minisketch::ImplementationSupported(bits, impl)) return;

    Minisketch sketch_a{Minisketch(bits, impl, capacity)};
    Minisketch sketch_b{Minisketch(bits, impl, capacity)};
    FUZZ_CHECK(sketch_a && sketch_b);
    sketch_a.SetSeed(fuzzed_data_provider.ConsumeIntegral<uint64_t>());
    sketch_b.SetSeed(fuzzed_data_provider.ConsumeIntegral<uint64_t>());

    // Fill two sets and keep the difference in a map.
    std::map<uint64_t, bool> diff;
    for (int i = 0; i < 10000 && fuzzed_data_provider.ConsumeBool(); ++i) {
        const auto entry{fuzzed_data_provider.ConsumeIntegralInRange<uint64_t>(1, max_entry)};
        const auto keep_diff{[&] {
            bool& mut{diff[entry]};
            mut = !mut;
        }};
        switch (fuzzed_data_provider.ConsumeIntegralInRange<int>(0, 2)) {
        case 0:
            sketch_a.Add(entry);
            keep_diff();
            break;
        case 1:
            sketch_b.Add(entry);
            keep_diff();
            break;
        case 2:
            sketch_a.Add(entry);
            sketch_b.Add(entry);
            break;
        }
    }
    const auto num_diff{std::accumulate(diff.begin(), diff.end(), size_t{0}, [](auto n, const auto& e) { return n + e.second; })};

    Minisketch sketch_ar{Minisketch(bits, impl, capacity)};
    Minisketch sketch_br{Minisketch(bits, impl, capacity)};
    sketch_ar.SetSeed(fuzzed_data_provider.ConsumeIntegral<uint64_t>());
    sketch_br.SetSeed(fuzzed_data_provider.ConsumeIntegral<uint64_t>());

    sketch_ar.Deserialize(sketch_a.Serialize());
    sketch_br.Deserialize(sketch_b.Serialize());

    Minisketch sketch_diff{std::move(fuzzed_data_provider.ConsumeBool() ? sketch_a : sketch_ar)};
    sketch_diff.Merge(fuzzed_data_provider.ConsumeBool() ? sketch_b : sketch_br);

    if (capacity >= num_diff) {
        const auto max_elements{fuzzed_data_provider.ConsumeIntegralInRange<size_t>(num_diff, capacity)};
        const auto dec{sketch_diff.Decode(max_elements)};
        FUZZ_CHECK(dec.has_value());
        FUZZ_CHECK(dec->size() == num_diff);
        for (auto d : *dec) {
            FUZZ_CHECK(diff.at(d));
        }
    }
}
