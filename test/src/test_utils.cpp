#include "test_utils.hpp"

namespace pef_test_utils
{
// helper function for NextGEQ test
uint64_t scan_nextGEQ2(sdsl::bit_vector &bv, uint64_t x) {
    uint64_t res = 0;
    for (uint64_t i = x; i < bv.size(); i++) {
        if (bv[i] == 1) {
            res = i;
            break;
        }
    }
    return res;
}
}
