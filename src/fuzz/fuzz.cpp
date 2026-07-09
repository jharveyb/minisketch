/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include "fuzz.h"

#include <cstring>
#include <map>
#include <string>

namespace {

std::map<std::string, FuzzTargetFn>& Targets() {
    static std::map<std::string, FuzzTargetFn> targets;
    return targets;
}

FuzzTargetFn g_target = nullptr;

} // namespace

bool RegisterFuzzTarget(const char* name, FuzzTargetFn fn) {
    Targets().emplace(name, fn);
    return true;
}

extern "C" int LLVMFuzzerInitialize(int*, char***) {
    const char* name = getenv("FUZZ");
    if (name != nullptr) {
        auto it = Targets().find(name);
        if (it != Targets().end()) g_target = it->second;
    }
    if (g_target == nullptr) {
        fprintf(stderr, "Select a fuzz target with the FUZZ environment variable, e.g. FUZZ=decode %s\n",
                name == nullptr ? "" : "(unknown target)");
        fprintf(stderr, "Available targets:\n");
        for (const auto& entry : Targets()) fprintf(stderr, "  %s\n", entry.first.c_str());
        exit(1);
    }
    return 0;
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    g_target(data, size);
    return 0;
}
