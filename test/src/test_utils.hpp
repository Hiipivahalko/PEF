#ifndef _TEST_UTILS_HPP
#define _TEST_UTILS_HPP

#include <sdsl/bit_vectors.hpp>

namespace pef_test_utils
{
// helper function for NextGEQ test
uint64_t scan_nextGEQ2(sdsl::bit_vector &bv, uint64_t x);
}
#endif
