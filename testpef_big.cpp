#include "pef_vector_unif.hpp"
#include "pef_vector_opt.hpp"
#include "fstream"
#include "iostream"

#define N 196433254820

int main() {
  std::vector<uint64_t> pb;
  bit_vector b(N);
  uint64_t pos;

  ifstream rf("/data/bitvectors/ii/gov2/url/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);

  if(!rf){
    std::cout << "Cannot open file!\n";
    return 1;
  }
  
  while(!rf.eof()){
    rf.read((char*) &pos, sizeof(uint64_t));
    pb.push_back(pos);
    b[pos] = 1;
  }
  rf.close();
  
  sd_vector<> sd(b);
  cout << "Size PEF SD_VECTOR " << size_in_bytes(sd) << " " << ((double) (size_in_bytes(sd) * 8)) / ((double)pb.size()) << "\n";
  pef_vector_unif<1024, rank_support_v5<1>, select_support_mcl<1>> pef_unif(pb, N);
  cout << "Size PEF UNIF BIT VECTOR POS ONES " << pef_unif.size_in_bytes() << " " << ((double) (pef_unif.size_in_bytes() * 8)) / ((double)pb.size()) << "\n";
  pef_vector_opt<rank_support_v5<1>, select_support_mcl<1>> pef_opt(pb, N, 0.03, 0.3);
  cout << "Size PEF OPT BIT VECTOR POS ONES " << pef_opt.size_in_bytes() << " " << ((double) (pef_opt.size_in_bytes() * 8)) / ((double)pb.size()) << "\n";
  return 0;
}