#include <iostream>
#include <string>
#include <sdsl/bit_vectors.hpp>
#include <vector>
#include <iomanip>
#include <chrono>
//#include "pef_vector_unif.hpp"
#include <pef/pef_vector_opt.hpp>
//#include "pef_vector_opt_vigna.hpp"



int main(int argc, char *argv[]) {

	if (argc < 2) {
        std::cout << ">> Program usage: " << argv[0] << " <input_file> \n";
        exit(1);
    }

    std::cerr << ">> Program START time: " << std::endl;
    uint64_t postings_list_size = 0;
    std::cerr << ">> loading input bitvector...";
    sdsl::bit_vector bv;
    sdsl::load_from_file(bv, argv[1]);
    bv[bv.size()-1] = 1;
    std::cerr << "DONE\n";
    //postings_list_size = sdsl::rank_support_v5<1>(&bv)(bv.size());


    // creating test query values
    uint64_t n = bv.size() > 1000000 ? 1000000 : bv.size();
    uint64_t u = bv.size()-2;
    std::vector<uint64_t> benchmark_quesries(n);
    srand(0);
    for (uint32_t i = 0; i < n; i++) benchmark_quesries[i] = rand() % u + 1;


    std::cerr << ">> building pef...";
    pef_vector_opt<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef(bv, 0.03, 0.3);
    //pef_vector_opt<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef(bv, 0.5, 0.5);
    std::cerr << "done\n";

    std::cerr << "testing pef nextGEQ -> ...";
    for (uint64_t x = 1; x < bv.size(); x++) {
        //std::cerr << "x:" << x <<"\n";
        uint64_t succ = pef.nextGEQ(x);
    }
    std::chrono::duration<double, std::micro> ms_double;
        double time_avg = 0;
    double time_avg_total = 0;
    time_avg = 0;
    time_avg_total = 0;
    for (uint32_t j = 0; j < 3; j++) {
      time_avg = 0;
      for (uint64_t i = 0; i < benchmark_quesries.size(); i++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        pef.nextGEQ(benchmark_quesries[i]);
        auto t2 = std::chrono::high_resolution_clock::now();
        ms_double = t2 - t1;
        time_avg += ms_double.count();
      }
      time_avg_total += (time_avg / benchmark_quesries.size());
    }
    std::cerr << "DONE\n";
    std::cout << "avg nextGEQ:" << (time_avg_total / 3) << " mirco sec\n";

    std::cout << "pef size:" << ((double)pef.size_in_bits_nextGEQ() / 5055078461) << "\n";


    std::cerr << ">> Program END" << std::endl;

    return 0;
}
