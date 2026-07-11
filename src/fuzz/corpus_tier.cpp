/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Corpus tiering helper for the decode/roundtrip/reconcile fuzz targets.
 *
 * Decode cost is dominated by the `capacity` parameter (quadratic in it), so a
 * corpus can be split into a fast tier (small capacity, always loaded) and a
 * slow tier (high capacity, loaded on demand). This tool prints the capacity
 * each corpus file decodes to, mirroring the target's FuzzedDataProvider
 * consumption so routing stays exactly in sync with the harness. Policy (the
 * threshold and where files go) lives in tools/update_seed_corpus.sh.
 *
 *   corpus-tier <decode|roundtrip|reconcile> <corpus-dir>
 *
 * prints one "<capacity>\t<filename>" line per regular file in the directory.
 */

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#include "FuzzedDataProvider.h"

namespace fs = std::filesystem;

namespace {

/** The capacity a target derives from an input, matching decode.cpp /
 * roundtrip.cpp (bits via ConsumeIntegralInRange) and reconcile.cpp (bits via
 * ConsumeBool). The subsequent log-uniform capacity draw is identical across
 * all three. */
size_t CapacityOf(const uint8_t* data, size_t size, const std::string& target) {
    FuzzedDataProvider provider(data, size);
    if (target == "reconcile") {
        provider.ConsumeBool();
    } else {
        provider.ConsumeIntegralInRange<uint32_t>(2, 64);
    }
    size_t cap_lo = size_t{1} << provider.ConsumeIntegralInRange<int>(0, 10);
    return provider.ConsumeIntegralInRange<size_t>(cap_lo, std::min<size_t>(2 * cap_lo - 1, 1024));
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s <decode|roundtrip|reconcile> <corpus-dir>\n", argv[0]);
        return 2;
    }
    const std::string target = argv[1];
    if (target != "decode" && target != "roundtrip" && target != "reconcile") {
        fprintf(stderr, "%s: unknown target '%s'\n", argv[0], target.c_str());
        return 2;
    }
    std::vector<uint8_t> buf(1 << 20);
    for (const auto& entry : fs::directory_iterator(argv[2])) {
        if (!entry.is_regular_file()) continue;
        FILE* f = fopen(entry.path().c_str(), "rb");
        if (!f) continue;
        size_t n = fread(buf.data(), 1, buf.size(), f);
        fclose(f);
        printf("%zu\t%s\n", CapacityOf(buf.data(), n, target), entry.path().filename().string().c_str());
    }
    return 0;
}
