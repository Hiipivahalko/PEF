#include "pef_vector_unif.hpp"
#include "pef_vector_opt.hpp"
#include "pef_vector_opt_vigna.hpp"

#define N 15000000
//#define N 14


int main() {
  sdsl::bit_vector b(N);// = {0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1};
  //bit_vector b = {0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,1,1,1};
  //bit_vector b =  {0,0,1,1,0,1,0,1,0,1,0,1,0,0};

  //bit_vector test = {1,0,1,1,0,1,0,1,0,1,0,1,0,0};
  sdsl::util::set_random_bits(b);
  //sdsl::util::set_to_value(b, 1);

  //for (int i = 0; i < N; ++i){
  //  if(i % 5 == 0 || i % 5 == 1){
  //    b[i] = 1;
  //  }
  //}
  sdsl::sd_vector<> sd(b);
  //sd_vector<> sd_t(test);
  sdsl::rank_support_sd<1> rank_sd(&sd);    
  //rank_support_sd<1> rank_test(&sd_t);  
  //cout << rank_test(2) << "\n";  
  sdsl::select_support_sd<1> select_sd(&sd); 

  uint64_t n = rank_sd(N);
  cout << "---------------- PEF UNIF TESTS ----------------\n";
    pef_vector_unif<100, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef_unif(b);

    cout << "# Select\n";
    bool pass_tests = 1;
    for (uint64_t i = 1; i <= n; i++) {
      if (pef_unif.select(i) != select_sd(i)){
        pass_tests = 0;
        cout << pef_unif.select(i) << " " << select_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";

    cout << "# Rank\n";
    for (uint64_t i = 0; i <= N; i++) {
      if (pef_unif.rank(i) != rank_sd(i)){
        pass_tests = 0;
        cout << pef_unif.rank(i) << " " << rank_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";

    pass_tests = 1;

    cout << "# Access\n";
    pass_tests = 1;
    for (uint64_t i = 0; i < N; i++) {
      //cout << "bit " << i << " " << b[i] << "\n";
      if (pef_unif[i] != b[i]){
        pass_tests = 0;
        cout << pef_unif[i] << " " << b[i] << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
 cout << "-------- PEF UNIF TESTS POSITION OF ONES --------\n";
    std::vector<uint64_t> pb;

    for(uint64_t i = 0; i < N; i++){
      if(b[i] == 1) pb.push_back(i);
    }
    pef_vector_unif<100, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef_unif2(pb, N);
    cout << "# Select\n";
    pass_tests = 1;
    for (uint64_t i = 1; i <= n; i++) {
      if (pef_unif2.select(i) != select_sd(i)){
        pass_tests = 0;
        cout << pef_unif2.select(i) << " " << select_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
    pass_tests = 1;
    cout << "# Rank\n";
    for (uint64_t i = 0; i <= N; i++) {
      if (pef_unif2.rank(i) != rank_sd(i)){
        pass_tests = 0;
        cout << pef_unif2.rank(i) << " " << rank_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
    pass_tests = 1;

    cout << "# Access\n";
    pass_tests = 1;
    for (uint64_t i = 0; i < N; i++) {
      if (pef_unif2[i] != b[i]){
        pass_tests = 0;
        cout << pef_unif2[i] << " " << b[i] << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
 cout << "---------------- PEF OPTI TESTS ----------------\n";

    pef_vector_opt<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef_opt(b, 0.03, 0.3);
    cout << "# Select\n";
    pass_tests = 1;
    for (uint64_t i = 1; i <= n; i++) {
      if (pef_opt.select(i) != select_sd(i)){
        pass_tests = 0;
        cout << pef_opt.select(i) << " " << select_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
    pass_tests = 1;
    cout << "# Rank\n";
    for (uint64_t i = 0; i <= N; i++) {
      if (pef_opt.rank(i) != rank_sd(i)){
        pass_tests = 0;
        cout << pef_opt.rank(i) << " " << rank_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
    pass_tests = 1;

    cout << "# Access\n";
    pass_tests = 1;
    for (uint64_t i = 0; i < N; i++) {
      if (pef_opt[i] != b[i]){
        pass_tests = 0;
        cout << pef_opt[i] << " " << b[i] << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";

  cout << "-------- PEF OPTI TESTS POSITION OF ONES --------\n";

    pef_vector_opt<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef_opt2(pb, N, 0.03, 0.3);
    cout << "# Select\n";
    pass_tests = 1;
    for (uint64_t i = 1; i <= n; i++) {
      if (pef_opt2.select(i) != select_sd(i)){
        pass_tests = 0;
        cout << pef_opt2.select(i) << " " << select_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
    pass_tests = 1;
    cout << "# Rank\n";
    for (uint64_t i = 0; i <= N; i++) {
      if (pef_opt2.rank(i) != rank_sd(i)){
        pass_tests = 0;
        cout << pef_opt2.rank(i) << " " << rank_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
    pass_tests = 1;

    cout << "# Access\n";
    pass_tests = 1;
    for (uint64_t i = 0; i < N; i++) {
      if (pef_opt2[i] != b[i]){
        pass_tests = 0;
        cout << pef_opt2[i] << " " << b[i] << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
  cout << "-------- PEF OPTI TESTS POSITION OF ONES --------\n";

    pef_vector_opt_vigna<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef_opt_vigna_2(pb, N, 0.03, 0.3);
    cout << "# Select\n";
    pass_tests = 1;
    for (uint64_t i = 1; i <= n; i++) {
      if (pef_opt_vigna_2.select(i) != select_sd(i)){
        pass_tests = 0;
        cout << pef_opt_vigna_2.select(i) << " " << select_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
    pass_tests = 1;
    cout << "# Rank\n";
    for (uint64_t i = 0; i <= N; i++) {
      if (pef_opt_vigna_2.rank(i) != rank_sd(i)){
        pass_tests = 0;
        cout << pef_opt_vigna_2.rank(i) << " " << rank_sd(i) << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
    pass_tests = 1;

    cout << "# Access\n";
    pass_tests = 1;
    for (uint64_t i = 0; i < N; i++) {
      if (pef_opt_vigna_2[i] != b[i]){
        pass_tests = 0;
        cout << pef_opt_vigna_2[i] << " " << b[i] << "\n";
      }
    }
    if(pass_tests) cout << "No difference!\n";
  cout << "------------------------------------------------\n";

  cout << "Size EF " << ((double)size_in_bytes(sd) * 8) / pb.size() << "\n";
  cout << "Size PEF UNIF BIT VECTOR " << ((double)pef_unif.size_in_bytes() * 8) / pb.size() << "\n";
  cout << "Size PEF UNIF POS VECTOR " << ((double)pef_unif2.size_in_bytes() * 8) / pb.size() << "\n";
  cout << "Size PEF OPT BIT VECTOR " << ((double)pef_opt.size_in_bytes() * 8) / pb.size() << "\n";
  cout << "Size PEF OPT POS VECTOR " << ((double)pef_opt2.size_in_bytes() * 8) / pb.size() << "\n";
  cout << "Size PEF OPT VIGNA POS VECTOR " << ((double)pef_opt_vigna_2.size_in_bytes() * 8) / pb.size() << "\n";

  return 0;
}

