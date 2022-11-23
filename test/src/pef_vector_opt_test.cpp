#include <gtest/gtest.h>
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <tuple>
#include <math.h>
#include <string>

#include "test_utils.hpp"

#include <pef/pef_vector_opt.hpp>


TEST(PEFNextGEQ, smallTestSample1) {
    sdsl::bit_vector bv;
    sdsl::load_from_file(bv, "../data/small_sample_as_bv.dat");
    bv[bv.size()-1] = 1;
    std::cerr << ">> building pef...";
    //pef_vector_opt<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef(bv, 0.5, 0.5);
    pef_vector_opt<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef(bv, 0.03, 0.3);
    std::cerr << "DONE\n";

    for (int i = 1; i < bv.size(); i++) {
        uint64_t res_scan = pef_test_utils::scan_nextGEQ2(bv, i);
        //std::cerr << "res_Scan:" << res_scan << "\n";
        ASSERT_EQ(pef.nextGEQ(i), res_scan);
    }
}

TEST(PEFNextGEQ, smallTestSample2) {
    sdsl::bit_vector bv;
    sdsl::load_from_file(bv, "../data/small_sample_as_bv.dat");
    bv[bv.size()-1] = 1;
    std::cerr << ">> building pef...";
    pef_vector_opt<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef(bv, 0.02, 0.2);
    std::cerr << "DONE\n";

    for (int i = 1; i < bv.size(); i++) {
        uint64_t res_scan = pef_test_utils::scan_nextGEQ2(bv, i);
        //std::cerr << "res_Scan:" << res_scan << "\n";
        ASSERT_EQ(pef.nextGEQ(i), res_scan);
    }
}
