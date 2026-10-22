/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** oracle_dump's clmul field instantiations (a separate TU from
 * oracle_dump.cpp: the generic and clmul headers both define `Field` in the
 * anonymous namespace). */

#include <string>

#ifdef HAVE_CLMUL
#include "fields/clmul_1byte.cpp"
#include "fields/clmul_2bytes.cpp"
#include "fields/clmul_4bytes.cpp"
#endif

#include "oracle_dump_impl.h"

namespace oracle_dump {

std::string RunClmul(int impl, int bits, const Case& c) {
#ifdef HAVE_CLMUL
#define ORACLE_CASE(B, T) case B: return Run<T>(impl, B, c);
    if (impl == 1) {
        switch (bits) {
#ifdef ENABLE_FIELD_INT_5
        ORACLE_CASE(5, Field5)
#endif
#ifdef ENABLE_FIELD_INT_8
        ORACLE_CASE(8, Field8)
#endif
#ifdef ENABLE_FIELD_INT_10
        ORACLE_CASE(10, Field10)
#endif
#ifdef ENABLE_FIELD_INT_11
        ORACLE_CASE(11, Field11)
#endif
#ifdef ENABLE_FIELD_INT_12
        ORACLE_CASE(12, Field12)
#endif
#ifdef ENABLE_FIELD_INT_13
        ORACLE_CASE(13, Field13)
#endif
#ifdef ENABLE_FIELD_INT_14
        ORACLE_CASE(14, Field14)
#endif
#ifdef ENABLE_FIELD_INT_16
        ORACLE_CASE(16, Field16)
#endif
#ifdef ENABLE_FIELD_INT_32
        ORACLE_CASE(32, Field32)
#endif
        }
    } else if (impl == 2) {
        switch (bits) {
#ifdef ENABLE_FIELD_INT_2
        ORACLE_CASE(2, FieldTri2)
#endif
#ifdef ENABLE_FIELD_INT_3
        ORACLE_CASE(3, FieldTri3)
#endif
#ifdef ENABLE_FIELD_INT_4
        ORACLE_CASE(4, FieldTri4)
#endif
#ifdef ENABLE_FIELD_INT_5
        ORACLE_CASE(5, FieldTri5)
#endif
#ifdef ENABLE_FIELD_INT_6
        ORACLE_CASE(6, FieldTri6)
#endif
#ifdef ENABLE_FIELD_INT_7
        ORACLE_CASE(7, FieldTri7)
#endif
#ifdef ENABLE_FIELD_INT_9
        ORACLE_CASE(9, FieldTri9)
#endif
#ifdef ENABLE_FIELD_INT_10
        ORACLE_CASE(10, FieldTri10)
#endif
#ifdef ENABLE_FIELD_INT_11
        ORACLE_CASE(11, FieldTri11)
#endif
#ifdef ENABLE_FIELD_INT_12
        ORACLE_CASE(12, FieldTri12)
#endif
#ifdef ENABLE_FIELD_INT_14
        ORACLE_CASE(14, FieldTri14)
#endif
#ifdef ENABLE_FIELD_INT_15
        ORACLE_CASE(15, FieldTri15)
#endif
        }
    }
#undef ORACLE_CASE
#endif
    (void)impl; (void)bits; (void)c;
    return "";
}

} // namespace oracle_dump
