#include "sux/sux/bits/EliasFano.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#define N 10

inline int lambda_safe(uint64_t word) { return word == 0 ? -1 : 63 ^ __builtin_clzll(word); }

int main(){
  sdsl::bit_vector b(N);// = {0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1};
  sdsl::util::set_random_bits(b);

  std::vector< uint64_t > pb;
  for(uint64_t i = 0; i < N; i++){
    if(b[i]) pb.push_back(i); 
  }
  sux::bits::EliasFano<> ef(pb, N);
  std::cout << ef.select(1) << "\n";
  std::cout << ef.bitCount() << "\n";  
  std::cout << ef.bitCount() / 8 << "\n";  

  uint64_t nums_ones = pb.size();
  uint64_t nums_bits = N;
  uint64_t L = (nums_ones == 0 ? 0 : std::max(lambda_safe(nums_bits / nums_ones), 0));
  uint64_t size_lower = ((nums_ones * L + 63) / 64 + 2 * (L == 0)) * sizeof(uint64_t);
  uint64_t upper_bits = nums_ones + (nums_bits >> L);
  uint64_t size_upper = ((upper_bits + 1 + 63) / 64) * sizeof(uint64_t);

  uint64_t log2_ones_per_inventory = 10;
  uint64_t ones_per_inventory = 1 << 10;
  uint64_t log2_longwords_per_subinventory = 2;
  uint64_t longwords_per_subinventory = 1 << log2_longwords_per_subinventory;
  uint64_t size_select_upper = (((nums_ones + ones_per_inventory - 1) / ones_per_inventory) * (longwords_per_subinventory + 1) + 1) * sizeof(uint64_t);
 
  std::cout << upper_bits << " " << nums_ones << "\n"; 
  uint64_t nums_zeros = upper_bits - nums_ones;
  std::cout << nums_zeros << "\n";
  uint64_t log2_zeros_per_inventory = 10;
  uint64_t zeros_per_inventory = 1 << log2_zeros_per_inventory;
  uint64_t size_selectz_upper = (((nums_zeros + zeros_per_inventory - 1) / zeros_per_inventory) * (longwords_per_subinventory + 1) + 1) * sizeof(uint64_t);

  std::cout << size_upper << " ";
  std::cout << size_lower << " ";
  std::cout << size_select_upper << " ";
  std::cout << size_selectz_upper << " ";

  std::cout << sizeof(sux::bits::EliasFano<>) + size_lower + size_upper + size_select_upper + size_selectz_upper << "\n";
  std::cout << 8*(sizeof(sux::bits::EliasFano<>) + size_lower + size_upper + size_select_upper + size_selectz_upper) << "\n";
 
  sdsl::sd_vector<> sd(pb.begin(), pb.end());
  std::cout << sdsl::size_in_bytes(sd) << "\n";
  std::cout << sdsl::size_in_bytes(b) << "\n"; 
  return 0;
}
