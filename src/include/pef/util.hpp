#ifndef PEF_UTILS
#define PEF_UTILS

#include <vector>
#include <utility>
//#include <sux/bits/EliasFano.hpp>

using namespace std;

class cost_window {
  public:
    std::vector<uint64_t>::iterator _start_it;
    std::vector<uint64_t>::iterator _end_it;
    uint64_t _start = 0;
    uint64_t _end = 0; // end-th position is not in the current window
    uint64_t _min_p = 0; // element that preceed the first element of the window
    uint64_t _max_p = 0;

    uint64_t _cost_upper_bound; // The maximum cost for this window

    cost_window(std::vector<uint64_t>::iterator begin, uint64_t __cost_upper_bound) {
      _start_it = begin;
      _end_it = begin;
      _min_p = *begin;
      _max_p = 0;
      _cost_upper_bound = __cost_upper_bound;
    }

    uint64_t universe() const {
      return _max_p - _min_p + 1;
    }

    uint64_t size() const {
      return _end - _start;
    }

    void advance_start() {
      _min_p = *_start_it + 1;
      ++_start;
      ++_start_it;
    }

    void advance_end() {
      _max_p = *_end_it;
      ++_end;
      ++_end_it;
    }

    uint64_t start() {
      return _start;
    }

    uint64_t end() {
      return _end;
    }

    uint64_t cost_upper_bound() {
      return _cost_upper_bound;
    }

    void print(uint64_t i) {
      cout << "window " << i << ": start "
           << _start << " end "
           << _end << " min_p "
           << _min_p << " max_p "
           << _max_p << " u "
           << universe() << " n "
           << size() << " cost "
           << _cost_upper_bound << "\n";
    }
};

static const uint64_t type_bits = 1; // all_ones is implicit

// from succint repository
static const uint8_t debruijn64_mapping[64] = {
  63,  0, 58,  1, 59, 47, 53,  2,
  60, 39, 48, 27, 54, 33, 42,  3,
  61, 51, 37, 40, 49, 18, 28, 20,
  55, 30, 34, 11, 43, 14, 22,  4,
  62, 57, 46, 52, 38, 26, 32, 41,
  50, 36, 17, 19, 29, 10, 13, 21,
  56, 45, 25, 31, 35, 16,  9, 12,
  44, 24, 15,  8, 23,  7,  6,  5
};

static const uint64_t debruijn64 = 0x07EDD5E59A4E28C2ULL;

// return the position of the single bit set in the word x
inline uint8_t bit_position(uint64_t x){
    assert(__builtin_popcountll(x) == 1);
    return debruijn64_mapping[(x * debruijn64) >> 58];
}

inline uint8_t msb(uint64_t x, unsigned long& ret){
  if (!x) {
    return false;
  }

  // right-saturate the word
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;

  // isolate the MSB
  x ^= x >> 1;
  ret = bit_position(x);

  return true;
}

inline uint8_t msb(uint64_t x){
  unsigned long ret = -1U;
  msb(x, ret);
  return (uint8_t)ret;
}



inline uint64_t ceil_log2(uint64_t x) {
  return (x > 1) ? msb(x - 1) + 1 : 0;
}

// bitsize calculator for plain bit vector in https://github.com/ot/partitioned_elias_fano
uint64_t bitsize_plain_bitvector(uint64_t universe, uint64_t n){
  // global parameters
  uint64_t rb_log_rank1_sampling = 9; // cuidado con esto
  uint64_t rb_log_sampling1 = 8; // cuidado con esto

  uint64_t log_rank1_sampling = rb_log_rank1_sampling; // cuidado con esto
  uint64_t log_sampling1 = rb_log_sampling1; // cuidado con esto

  uint64_t rank1_sample_size = ceil_log2(n + 1);
  uint64_t pointer_size = ceil_log2(universe);
  uint64_t rank1_samples = universe >> log_rank1_sampling;
  uint64_t pointers1 = n >> log_sampling1;
  uint64_t rank1_samples_offset = 0;
  uint64_t pointers1_offset = rank1_samples_offset + rank1_samples * rank1_sample_size;
  uint64_t bits_offset = pointers1_offset + pointers1 * pointer_size;
  return bits_offset + universe;
}

// bitsize calculator for elias fano in https://github.com/ot/partitioned_elias_fano
uint64_t bitsize_elias_fano(uint64_t universe, uint64_t n){
  //global parameters
  uint64_t ef_log_sampling0 = 9; // cuidado con esto
  uint64_t ef_log_sampling1 = 8; // cuidado con esto
  uint64_t log_sampling0 = ef_log_sampling0;
  uint64_t log_sampling1 = ef_log_sampling1;
  uint64_t lower_bits = universe > n ? msb(universe / n) : 0;
  uint64_t mask = (((uint64_t) 1) << lower_bits) - 1;
  uint64_t higher_bits_length = n + (universe >> lower_bits) + 2;
  uint64_t pointer_size = ceil_log2(higher_bits_length);
  uint64_t pointers0 = (higher_bits_length - n) >> log_sampling0;
  uint64_t pointers1 = (n >> log_sampling1);
  uint64_t pointers0_offset = 0;
  uint64_t pointers1_offset = (pointers0_offset + pointers0 * pointer_size);
  uint64_t higher_bits_offset = pointers1_offset + pointers1 * pointer_size;
  uint64_t lower_bits_offset = higher_bits_offset + higher_bits_length;

  return lower_bits_offset + n * lower_bits;
}

inline int lambda_safe(uint64_t word) { return word == 0 ? -1 : 63 ^ __builtin_clzll(word); }

uint64_t bitsize_ef_vigna(uint64_t universe, uint64_t n) {
  uint64_t nums_ones = n;
  uint64_t nums_bits = universe;

  // size in bytes of lower attribute in EliasFano class of Vigna
  uint64_t L = (nums_ones == 0 ? 0 : std::max(lambda_safe(nums_bits / nums_ones), 0));
  uint64_t size_lower = ((nums_ones * L + 63) / 64 + 2 * (L == 0)) * sizeof(uint64_t);

  // size in bytes of upper attribute in EliasFano class of Vigna
  uint64_t upper_bits = nums_ones + (nums_bits >> L);
  uint64_t size_upper = ((upper_bits + 1 + 63) / 64) * sizeof(uint64_t);

  // this two variables are usefull to know the size of SimpleSelectHalf and SimpleSelectZeroHalf
  // two attributes in EliasFano class of vigna
  uint64_t log2_longwords_per_subinventory = 2;
  uint64_t longwords_per_subinventory = 1 << log2_longwords_per_subinventory;

  // size in bytes of select_upper attribute in EliasFano class of Vigna
  uint64_t log2_ones_per_inventory = 10;
  uint64_t ones_per_inventory = 1 << 10;
  uint64_t size_select_upper = 0;//(((nums_ones + ones_per_inventory - 1) / ones_per_inventory) * (longwords_per_subinventory + 1) + 1) * sizeof(uint64_t);

  // size in bytes of selectz_upper attribute in EliasFano class of Vigna
  uint64_t nums_zeros = upper_bits - nums_ones;
  uint64_t log2_zeros_per_inventory = 10;
  uint64_t zeros_per_inventory = 1 << log2_zeros_per_inventory;
  uint64_t size_selectz_upper = 0;//(((nums_zeros + zeros_per_inventory - 1) / zeros_per_inventory) * (longwords_per_subinventory + 1) + 1) * sizeof(uint64_t);

  return 8 * (size_lower + size_upper + size_select_upper + size_selectz_upper);
}

void print(std::vector< uint64_t > &v){
  uint64_t i;
  for(i = 0; i < ((uint64_t) v.size()) - 1; i++){
    cout << v[i] << " ";
  }
  cout << v[i] << "\n";
}

uint64_t bitsize_selects_ef_vigna(uint64_t universe, uint64_t n){
  uint64_t nums_ones = n;
  uint64_t nums_bits = universe;

  // size in bytes of lower attribute in EliasFano class of Vigna
  uint64_t L = (nums_ones == 0 ? 0 : std::max(lambda_safe(nums_bits / nums_ones), 0));

  // size in bytes of upper attribute in EliasFano class of Vigna
  uint64_t upper_bits = nums_ones + (nums_bits >> L);

  // this two variables are usefull to know the size of SimpleSelectHalf and SimpleSelectZeroHalf
  // two attributes in EliasFano class of vigna
  uint64_t log2_longwords_per_subinventory = 2;
  uint64_t longwords_per_subinventory = 1 << log2_longwords_per_subinventory;

  // size in bytes of select_upper attribute in EliasFano class of Vigna
  uint64_t log2_ones_per_inventory = 10;
  uint64_t ones_per_inventory = 1 << 10;
  uint64_t size_select_upper = (((nums_ones + ones_per_inventory - 1) / ones_per_inventory) * (longwords_per_subinventory + 1) + 1) * sizeof(uint64_t);

  // size in bytes of selectz_upper attribute in EliasFano class of Vigna
  uint64_t nums_zeros = upper_bits - nums_ones;
  uint64_t log2_zeros_per_inventory = 10;
  uint64_t zeros_per_inventory = 1 << log2_zeros_per_inventory;
  uint64_t size_selectz_upper = (((nums_zeros + zeros_per_inventory - 1) / zeros_per_inventory) * (longwords_per_subinventory + 1) + 1) * sizeof(uint64_t);

  return 8 * (size_selectz_upper + size_select_upper);
}

uint64_t bitsize_final_ef_vigna(uint64_t u, uint64_t n){
  return bitsize_ef_vigna(u, n) + bitsize_selects_ef_vigna(u, n);
}

uint64_t bits_rank_support_v(uint64_t u, uint64_t n){
  uint64_t capacity = ((u + 63) >> 6) << 6;
  uint64_t size_int_vector = ((capacity >> 9) + 1);
  uint64_t ceil_log2_n = (n % 2 == 0 ? lambda_safe(n) : lambda_safe(n) + 1);
  return ceil_log2_n * size_int_vector + 64 * size_int_vector;
}

uint64_t bits_rank_support_v5(uint64_t u, uint64_t n){
  uint64_t capacity = ((u + 63) >> 6) << 6;
  uint64_t size_int_vector = ((capacity >> 11) + 1);
  uint64_t ceil_log2_n = (n % 2 == 0 ? lambda_safe(n) : lambda_safe(n) + 1);
  return ceil_log2_n * size_int_vector
         + 64 * size_int_vector;
}

uint64_t bits_select_support_mcl(uint64_t u, uint64_t n){
  uint64_t size_int_vector = (n + 4096 - 1) >> 12;
  uint64_t ceil_log2_u = (u % 2 == 0 ? lambda_safe(u) : lambda_safe(u) + 1);
  uint64_t m_superblock_bits = ceil_log2_u * size_int_vector;
  //uint64_t mini_or_long_bits = 64 * (size_int_vector / 64 + 1);

  return (sizeof(uint64_t)) * 8
         + m_superblock_bits;
         //+ mini_or_long_bits;
}

uint64_t bits_sd_vector(uint64_t u, uint64_t n){
  uint64_t ceil_log2_u = (u % 2 == 0 ? lambda_safe(u) : lambda_safe(u) + 1);
  uint64_t ceil_log2_n = (n % 2 == 0 ? lambda_safe(n) : lambda_safe(n) + 1);

  uint64_t bits_m_low = n * (ceil_log2_u - ceil_log2_n);

  uint64_t size_m_high = n + (1ULL << ceil_log2_n);
  uint64_t bits_m_high = 64 * (size_m_high / 64 + 1);

  uint64_t bits_m_high_1_select = bits_select_support_mcl(u, n);
  uint64_t bits_m_high_0_select = bits_select_support_mcl(u, u - n);

  return bits_m_low + bits_m_high
         + bits_m_high_1_select
         + bits_m_high_0_select;
}

//uint64_t bits_built_sd_vector(uint64_t u, uint64_t n, sdsl::sd_vector<> &sd){
//  uint64_t ceil_log2_u = (u % 2 == 0 ? lambda_safe(u) : lambda_safe(u) + 1);
//  uint64_t ceil_log2_n = (n % 2 == 0 ? lambda_safe(n) : lambda_safe(n) + 1);
//
//  uint64_t bits_m_low = n * (ceil_log2_u - ceil_log2_n);
//
//  uint64_t size_m_high = n + (1ULL << ceil_log2_n);
//  uint64_t bits_m_high = 64 * (size_m_high / 64 + 1);
//
//  uint64_t bits_m_high_1_select = bits_select_support_mcl(u, n);
//  uint64_t bits_m_high_0_select = bits_select_support_mcl(u, u - n);
//
//  return bits_m_low + bits_m_high
//         + bits_m_high_1_select
//         + bits_m_high_0_select
//         + sd.get_mini_or_low_bits_mcl();
//}
//
//uint64_t bits_built_sd_vector_no_select(uint64_t u, uint64_t n){
//  uint64_t ceil_log2_u = (u % 2 == 0 ? lambda_safe(u) : lambda_safe(u) + 1);
//  uint64_t ceil_log2_n = (n % 2 == 0 ? lambda_safe(n) : lambda_safe(n) + 1);
//
//  uint64_t bits_m_low = n * (ceil_log2_u - ceil_log2_n);
//
//  uint64_t size_m_high = n + (1ULL << ceil_log2_n);
//  uint64_t bits_m_high = 64 * (size_m_high / 64 + 1);
//
//  uint64_t bits_m_high_1_select = 0;//bits_select_support_mcl(u, n);
//  uint64_t bits_m_high_0_select = 0;//bits_select_support_mcl(u, u - n);
//
//  return bits_m_low + bits_m_high
//         + bits_m_high_1_select
//         + bits_m_high_0_select;
//         //+ sd.get_mini_or_low_bits_mcl();
//}

//uint64_t bits_built_select_support_mcl(uint64_t u, uint64_t n, sdsl::select_support_mcl<1> &select){
//  return bits_select_support_mcl(u, n) + select.get_mini_or_long_bits();
//}
//
uint64_t bits_bit_vector(uint64_t u, uint64_t n){
  uint64_t bits_bv = 64 * (u / 64 + 1);
  uint64_t bits_select = bits_select_support_mcl(u, n);
  uint64_t bits_rank = bits_rank_support_v(u, n);
  return bits_bv + bits_select + bits_rank;
}
//
//uint64_t bits_built_bit_vector(uint64_t u, uint64_t n, sdsl::select_support_mcl<1> &select){
//  uint64_t bits_bv = 64 * (u / 64 + 1);
//  uint64_t bits_select = bits_built_select_support_mcl(u, n, select);
//  uint64_t bits_rank = bits_rank_support_v(u, n);
//  return bits_bv + bits_select + bits_rank;
//}
//
//uint64_t bits_built_bit_vector_no_select(uint64_t u, uint64_t n){
//  uint64_t bits_bv = 64 * (u / 64 + 1);
//  uint64_t bits_select = 0;//bits_built_select_support_mcl(u, n, select);
//  uint64_t bits_rank = bits_rank_support_v(u, n);
//  return bits_bv + bits_select + bits_rank;
//}

// best bitsize function in https://github.com/ot/partitioned_elias_fano
uint64_t bitsize(uint64_t universe, uint64_t n){
  uint64_t best_cost;

  // the sequence has all 1s? --> 0 bits; otherwise, +infty
  best_cost  = (universe == n) ? 0 : uint64_t(-1);
  //cout << "universe: " << universe << " n: " << n << " ";

  uint64_t ef_cost = bitsize_elias_fano(universe, n) + type_bits;
  if (ef_cost < best_cost) {
    best_cost = ef_cost;
  }

  uint64_t rb_cost = bitsize_plain_bitvector(universe, n) + type_bits;
  if (rb_cost < best_cost) {
    best_cost = rb_cost;
  }
  //cout << ef_cost << " " << rb_cost << "\n";
  return best_cost;
}

uint64_t new_bitsize(uint64_t universe, uint64_t n){
  uint64_t best_cost;

  // the sequence has all 1s? --> 0 bits; otherwise, +infty
  best_cost  = (universe == n) ? 0 : uint64_t(-1);
  //cout << "universe: " << universe << " n: " << n << " ";

  uint64_t ef_cost = bits_sd_vector(universe, n) + type_bits;
  //std::cout << ef_cost << "\n";
  if (ef_cost < best_cost) {
    best_cost = ef_cost;
  }

  uint64_t rb_cost = bits_bit_vector(universe, n) + type_bits;
  //std::cout << rb_cost << "\n";
  if (rb_cost < best_cost) {
    best_cost = rb_cost;
  }
  //cout << ef_cost << " " << rb_cost << "\n";
  return best_cost;
}

uint64_t new_type_encoding(uint64_t universe, uint64_t n){
  uint64_t best_cost;
  uint64_t type;

  // the sequence has all 1s? --> 0 bits; otherwise, +infty
  best_cost  = (universe == n) ? 0 : uint64_t(-1);
  type = 0;
  //cout << "universe: " << universe << " n: " << n << " ";

  uint64_t ef_cost = bits_sd_vector(universe, n) + type_bits;
  //std::cout << ef_cost << "\n";
  if (ef_cost < best_cost) {
    best_cost = ef_cost;
    type = 1;
  }

  uint64_t rb_cost = bits_bit_vector(universe, n) + type_bits;
  //std::cout << rb_cost << "\n";
  if (rb_cost < best_cost) {
    best_cost = rb_cost;
    type = 2;
  }
  /*
  if(type == 1) {
    cout << "SD ";
  } else if(type == 2) {
    cout << "BIT ";
  }

  if(type != 0) cout << universe << " " << n << " " << best_cost << "\n";
*/
  //cout << ef_cost << " " << rb_cost << "\n";

  return type;
}

uint64_t type_encoding(uint64_t universe, uint64_t n){
  uint64_t best_cost;
  uint64_t type;

  // the sequence has all 1s? --> 0 bits; otherwise, +infty
  //if(universe == n + 1) cout << "universe: " << universe << " n: " << n << "\n";
  //if(universe == n + 2) cout << "universe: " << universe << " n: " << n << "\n";
  best_cost  = (universe == n) ? 0 : uint64_t(-1);
  type = 0;

  uint64_t ef_cost = bitsize_elias_fano(universe, n) + type_bits;
  if (ef_cost < best_cost) {
    best_cost = ef_cost;
    type = 1;
  }

  uint64_t rb_cost = bitsize_plain_bitvector(universe, n) + type_bits;
  if (rb_cost < best_cost) {
    best_cost = rb_cost;
    type = 2;
  }
  /*
  if(type == 1) {
    cout << "SD ";
  } else if(type == 2) {
    cout << "BIT ";
  }

  if(type != 0) cout << universe << " " << n << " " << best_cost << "\n";
  */
  return type;
}

// cost fun in github https://github.com/ot/partitioned_elias_fano
uint64_t fixed_cost = 192;

uint64_t cost_fun(uint64_t universe, uint64_t n) {
  // estimated best cost + fix cost
  return new_bitsize(universe, n) + fixed_cost; // config.fix_cost;  // esto último es 64? Ver, me parece que en general no.
};
#endif


