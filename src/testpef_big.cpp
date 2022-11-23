#include "pef_vector_unif.hpp"
#include "pef_vector_opt.hpp"
#include "fstream"
#include "iostream"

#define N 196433254820
#define BITS 80000000
#define GROUP_SIZE 20000000
int main() {
  std::vector<uint64_t> pb(GROUP_SIZE);
  //bit_vector b(N);
  uint64_t pos;
  ifstream rf("/data/bitvectors/ii/gov2/url/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);

  if(!rf){
    std::cout << "Cannot open file!\n";
    return 1;
  }

    //if(pb.size() > BITS) break;
  uint64_t i = 0;
  uint64_t last_one = 0;
  double sd_sum, pef_unif_sum, pef_opt_sum;
  sd_sum = pef_unif_sum = pef_opt_sum = 0;
  double times = 0;
  while(!rf.eof()){
    uint64_t readed = 0;
    uint64_t k = 0;
    while(!rf.eof()){
      rf.read((char*) &pos, sizeof(uint64_t));
      pb[k] = pos - last_one;
      k++;
      readed++;
      if(readed >= GROUP_SIZE) break;
    }
        
    sdsl::bit_vector b(pb[pb.size() - 1] + 1);
    sdsl::util::set_to_value(b, 0);
    for(uint64_t bit : pb) b[bit] = 1;
    
    sdsl::sd_vector<> sd(b);
    double size = ((double) (sdsl::size_in_bytes(sd) * 8)) / ((double)pb.size());
    cout << "Size PEF SD_VECTOR " << size << "\n";
    sd_sum += size;

    pef_vector_unif<4096, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef_unif(pb, pb[pb.size() - 1] + 1);
    size = ((double) (pef_unif.size_in_bytes() * 8)) / ((double)pb.size());
    cout << "Size PEF UNIF BIT VECTOR POS ONES " << size << "\n";
    pef_unif_sum += size;

    pef_vector_opt<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>, 2, 3> pef_opt(pb, pb[pb.size() - 1] + 1, 0.03, 0.3);
    size = ((double) (pef_opt.size_in_bytes() * 8)) / ((double)pb.size());
    cout << "Size PEF OPT BIT VECTOR POS ONES " << size << "\n";
    pef_opt_sum += size;
    cout << "\n";

    last_one = pb[pb.size() - 1];
    sdsl::util::clear(sd);
    sdsl::util::clear(b);

    times++;
    //if(times == 2)break;
  }
  rf.close();

  cout << "AVERAGE SIZE SD_VECTOR " << sd_sum / times << "\n"; 
  cout << "AVERAGE SIZE PEF_UNIF " << pef_unif_sum / times << "\n"; 
  cout << "AVERAGE SIZE PEF_OPT " << pef_opt_sum / times << "\n"; 
  
  return 0;
}
