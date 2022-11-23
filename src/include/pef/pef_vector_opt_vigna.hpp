#ifndef PEF_OPT_VIGNA
#define PEF_OPT_VIGNA

#include <vector>
#include <utility>
#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
//#include "sdsl-lite/include/sdsl/vectors.hpp"
//#include "sdsl-lite/include/sdsl/bit_vectors.hpp"
#include "sux/sux/bits/EliasFano.hpp"
#include <cinttypes>
#include <assert.h> 
#include "util.hpp"

using namespace std;

// best bitsize function in https://github.com/ot/partitioned_elias_fano
uint64_t bitsize_vigna(uint64_t universe, uint64_t n){
  uint64_t best_cost;

  // the sequence has all 1s? --> 0 bits; otherwise, +infty
  best_cost  = (universe == n) ? 0 : uint64_t(-1);
  //cout << "universe: " << universe << " n: " << n << " ";

  //sux::bits::EliasFano<> ef(ones, universe);
  //uint64_t ef_cost = ef.bitCount() + type_bits;
  uint64_t ef_cost = bitsize_ef_vigna(universe, n);
  if (ef_cost < best_cost) {
    best_cost = ef_cost;
  }

  //uint64_t rb_cost = bitsize_plain_bitvector(universe, n) + type_bits;
  uint64_t rb_cost = bits_bit_vector(universe, n) + type_bits;
  if (rb_cost < best_cost) {
    best_cost = rb_cost;
  }
  //cout << ef_cost << " " << rb_cost << "\n";
  return best_cost;
}

uint64_t type_encoding_vigna(uint64_t universe, uint64_t n){
  uint64_t best_cost;
  uint64_t type;

  // the sequence has all 1s? --> 0 bits; otherwise, +infty
  //if(universe == n + 1) cout << "universe: " << universe << " n: " << n << "\n"; 
  //if(universe == n + 2) cout << "universe: " << universe << " n: " << n << "\n"; 
  best_cost  = (universe == n) ? 0 : uint64_t(-1);
  type = 0;

  //sux::bits::EliasFano<> ef(ones, universe);
  //uint64_t ef_cost = ef.bitCount() + type_bits;
  uint64_t ef_cost = bitsize_ef_vigna(universe, n);
  if (ef_cost < best_cost) {
    best_cost = ef_cost;
    type = 1;
  }

  uint64_t rb_cost = bits_bit_vector(universe, n) + type_bits;
  if (rb_cost < best_cost) {
    best_cost = rb_cost;
    type = 2;
  }

  return type;
}

uint64_t cost_fun_vigna(uint64_t universe, uint64_t n) {
  // estimated best cost + fix cost
  return bitsize_vigna(universe, n) + fixed_cost; // config.fix_cost;  // esto último es 64? Ver, me parece que en general no.
};

template<class rank_support_c=sdsl::rank_support_scan<1>, class select_support_c=sdsl::select_support_scan<1>, uint64_t _fixed_cost=64 >
class pef_vector_opt_vigna {
  // last element of each block
  sux::bits::EliasFano<> L;

  // the accumulatives lengths of blocks
  sux::bits::EliasFano<> E;

  std::vector<void *> P;  // pointer to the block encoding
  sdsl::bit_vector B; // bit vector indicating the encoding of each block.
  // B[i] == 0 && P[i] == NULL -> All ones in block
  // B[i] == 0 && P[i] != NULL -> Plain bit vector encoding of the block
  // B[i] == 1                 -> Elias-Fano encoding ob the block
  std::vector<void *> block_select;
  std::vector<void *> block_rank;

  uint64_t u; // universe size  (bitvector length)
  uint64_t n; // set size (number of 1s in the bitvector)
  uint64_t nBlocks; 
    
  uint64_t size_in_bits_auxiliar;
  public:
    using size_type = size_t;

    pef_vector_opt_vigna() {;}

    ~pef_vector_opt_vigna() {;}

    uint64_t size(){
      return u;
    }

    uint64_t n_ones(){
      return n;
    }

    uint64_t size_in_bytes(){
      uint64_t u_L = L.size();
      uint64_t n_L = L.rank(u_L);
      uint64_t u_E = E.size();
      uint64_t n_E = E.rank(u_E);

      uint64_t only_encoding = (L.bitCount() - bitsize_selects_ef_vigna(u_L, n_L)) / 8  
                    + (E.bitCount() - bitsize_selects_ef_vigna(u_E, n_E)) / 8 
                    + sdsl::size_in_bytes(B) + 3 * sizeof(uint64_t) 
                    + nBlocks * sizeof(void *);
      uint64_t size = L.bitCount() / 8  
                    + E.bitCount() / 8 
                    + sdsl::size_in_bytes(B) + 3 * sizeof(uint64_t) 
                    + nBlocks * sizeof(void *);
      uint64_t ef_times = 0;
      uint64_t bit_times = 0;
      uint64_t all_ones_times = 0;
      for (uint64_t i=0; i < nBlocks; ++i) {
        if (B[i]) {
          size += (*(sux::bits::EliasFano<> *)P[i]).bitCount() / 8;
          uint64_t u_p = (*(sux::bits::EliasFano<> *)P[i]).size();
          uint64_t n_p = (*(sux::bits::EliasFano<> *)P[i]).rank(u_p);  
          only_encoding += ((*(sux::bits::EliasFano<> *)P[i]).bitCount() - bitsize_selects_ef_vigna(u_p, n_p)) / 8;
          if((*(sux::bits::EliasFano<> *)P[i]).bitCount() != bitsize_ef_vigna(u_p, n_p)) cout <<(*(sux::bits::EliasFano<> *)P[i]).bitCount() << " " << bitsize_ef_vigna(u_p, n_p) << "\n" ; 
          ef_times++;
        } else {
          if (P[i]) {
            only_encoding += sdsl::size_in_bytes(*(sdsl::bit_vector *)P[i]);
            size += sdsl::size_in_bytes(*(sdsl::bit_vector *)P[i])             
                  + sdsl::size_in_bytes(*(rank_support_c *)block_rank[i])
                  + sdsl::size_in_bytes(*(select_support_c *)block_select[i]);
            bit_times++;
	        } else {
            all_ones_times++;
          }
        }
      }
      cout << "EF " << ef_times << " BIT " << bit_times << " ALL_ONES " << all_ones_times << " " << (double) only_encoding / (double) n << " ";

      return size;
    }

    uint64_t size_in_bits_formula(){
      rank_support_c rank_B(&B);

      uint64_t u_L = L.size();
      uint64_t n_L = L.rank(u_L);
      
      uint64_t u_E = E.size();
      uint64_t n_E = E.rank(u_E);

      uint64_t n_B = rank_B(nBlocks);
      uint64_t u_B = nBlocks;

      uint64_t size = bitsize_final_ef_vigna(u_L, n_L) 
                    + bitsize_final_ef_vigna(u_E, n_E)
                    + ((u_B / 64 + 1) * 64);
                    + 3 * sizeof(uint64_t) * 8
                    + nBlocks * sizeof(void *) * 8;
      uint64_t size_L = bitsize_final_ef_vigna(u_L, n_L);
      uint64_t size_E = bitsize_final_ef_vigna(u_E, n_E);
      uint64_t size_B = ((u_B / 64 + 1) * 64);

      //cout << u << " " << n << " " 
      //     << size << " " << (double) size / n << " "
      //     << size_L << " " << (double) size_L / n << " "
      //     << size_E << " " << (double) size_E / n << " "
      //     << size_B << " " << (double) size_B / n << " ";
      
      uint64_t ef_times = 0;
      uint64_t ef_size = 0;
      
      uint64_t bit_times = 0;
      uint64_t bit_size = 0;
      
      uint64_t only_blocks = 0;

      uint64_t all_ones_times = 0;
      
      for (uint64_t i=0; i < nBlocks; ++i) {
        if (B[i]) {
          uint64_t u_p = (*(sux::bits::EliasFano<> *)P[i]).size();
          uint64_t n_p = (*(sux::bits::EliasFano<> *)P[i]).rank(u_p);
          size += bitsize_final_ef_vigna(u_p, n_p);
          ef_size += bitsize_final_ef_vigna(u_p, n_p);
          //only_blocks += bitsize_final_ef_vigna(u_p, n_p);
          only_blocks += bitsize_ef_vigna(u_p, n_p);

          //size += bitsize_ef_vigna(u_p, n_p);
          //ef_size += bitsize_ef_vigna(u_p, n_p);
          //only_blocks += bitsize_ef_vigna(u_p, n_p);

          //size += bitsize_elias_fano(u_p, n_p);
          //ef_size += bitsize_elias_fano(u_p, n_p);

	        ef_times++;
        } else {
          if (P[i]) {
	          uint64_t u_p = (*(sdsl::bit_vector *)P[i]).size();
            uint64_t n_p = (*(rank_support_c *)block_rank[i])(u_p);
            size += bits_built_bit_vector(u_p, n_p, *(select_support_c *)block_select[i]);
            bit_size += bits_built_bit_vector(u_p, n_p, *(select_support_c *)block_select[i]);
            //only_blocks += bits_built_bit_vector(u_p, n_p, *(select_support_c *)block_select[i]);
            only_blocks += bits_built_bit_vector_no_select(u_p, n_p);
            
            //size += bits_bit_vector(u_p, n_p);
            //bit_size += bits_bit_vector(u_p, n_p);
            //only_blocks += bits_bit_vector(u_p, n_p);
                       
            //size += bitsize_plain_bitvector(u_p, n_p);
            //bit_size += bitsize_plain_bitvector(u_p, n_p);
            bit_times++;
	        } else { 
            all_ones_times++;
          }
        }
      }
      //cout << "EF " << ef_times << " BIT " << bit_times << " ALL_ONES " << all_ones_times 
      //     << " " << ef_size << " " << bit_size << " "
      //     << (double) ef_size / n << " " << (double) bit_size / n << " "
      //     << only_blocks << " "
      //     << (double) only_blocks / n << "\n";
      return size;
    }
    
    size_t serialize(std::ostream& out, sdsl::structure_tree_node *v = nullptr, const std::string &name = "") const {
      return size_in_bits_auxiliar / 8;
    }

    std::pair<std::vector<uint64_t>, uint64_t> optimal_partition(std::vector<uint64_t> &ones_bv, double eps1, double eps2){
      //fixed_cost = 2 * ceil_log2(u) + ceil_log2(n);
      uint64_t single_block_cost = cost_fun_vigna(u, n);
      std::vector<uint64_t> min_cost(n+1, single_block_cost);
      min_cost[0] = 0;

      // create the required window: one for each power of approx_factor
      std::vector<cost_window> windows;
      uint64_t cost_lb = cost_fun_vigna(1, 1); // minimum cost

      uint64_t cost_bound = cost_lb;

      while (eps1 == 0 || cost_bound < cost_lb / eps1) {
        windows.emplace_back(ones_bv.begin(), cost_bound);
        if (cost_bound >= single_block_cost) break;
        cost_bound = cost_bound * (1 + eps2);
      }

      std::vector<uint64_t> path(n + 1, 0);
        
      for (uint64_t i = 0; i < n; i++) {
        size_t last_end = i + 1;
        for (auto& window: windows) {

          assert(window.start() == i);
          while (window._end < last_end) {
            window.advance_end();
          }

          uint64_t window_cost;
          while (true) {
            window_cost = cost_fun_vigna(window.universe(), window.size());
            if (min_cost[i] + window_cost < min_cost[window.end()]) {
              min_cost[window.end()] = min_cost[i] + window_cost;
              path[window.end()] = i;
            }
            last_end = window.end();
            if (window.end() == n) break;
            if (window_cost >= window.cost_upper_bound()) break;
            window.advance_end();
          }
          window.advance_start();
        }
      }

      uint64_t curr_pos = n;
      std::vector<uint64_t> partition;
      while (curr_pos != 0) {
        partition.push_back(curr_pos);
        curr_pos = path[curr_pos];
      }

      std::reverse(partition.begin(), partition.end());
      uint64_t cost_opt = min_cost[n];

      return {partition, cost_opt};
    }

    void add_block(sdsl::bit_vector &block_bv, int64_t type_encoding_block, uint64_t i){
      //all ones
      if(type_encoding_block == 0) {
        //std::cout << "* ALL ONES\n";
        B[i] = 0;
        P[i] = NULL; 
        block_select[i] = NULL;
        block_rank[i] = NULL;
      //elias-fano
      } else if(type_encoding_block == 1) {
        //std::cout << "* ELIAS FANO\n";
        B[i] = 1;
        std::vector<uint64_t> ones_bv;
        select_support_c select_1(&block_bv);
        rank_support_c rank_1(&block_bv);
        uint64_t u = block_bv.size();
        uint64_t n = rank_1(block_bv.size());
        //std::cout << "** SELECTING ONES FANO\n";
        //std::cout << u << " " << n << "\n";
        for(uint64_t i = 0; i < n; i++){
          ones_bv.push_back(select_1(i + 1));
        }
        //std::cout << "** OUT SELECTING ONES FANO\n";

        P[i] = new sux::bits::EliasFano<>(ones_bv, u);
        block_select[i] = NULL;
        block_rank[i] = NULL;
        ones_bv.clear();
        sdsl::util::clear(rank_1);
        sdsl::util::clear(select_1);
      //rb
      } else if(type_encoding_block == 2) {
        //std::cout << "* ADDING BIT VECTOR\n";
        B[i] = 0;
        P[i] = new sdsl::bit_vector(block_bv);
        block_select[i] = new select_support_c((sdsl::bit_vector *)P[i]); 
        block_rank[i] = new rank_support_c((sdsl::bit_vector *)P[i]); 
        //std::cout << "* FINISH ADDING BIT VECTOR\n";
      }
    }

    pef_vector_opt_vigna(sdsl::bit_vector &bv, double eps1, double eps2) {
      // setup fixed_cost 
      fixed_cost = _fixed_cost;

      select_support_c select_1(&bv);
      rank_support_c rank_1(&bv);
      sdsl::bit_vector block_bv;
      std::vector<uint64_t> elements_of_L;
      std::vector<uint64_t> elements_of_E;

      u = bv.size();
      n = rank_1(bv.size());
      std::vector<uint64_t> ones_bv;

      for (uint64_t i = 0; i < n; i++)
        ones_bv.push_back(select_1(i+1));
      
      //----- Optimal partition -----
      pair<std::vector<uint64_t>, uint64_t> ret = optimal_partition(ones_bv, eps1, eps2);
      std::vector<uint64_t> partition = ret.first;
      uint64_t cost_opt = ret.second;
      //-----------------------------

      nBlocks = partition.size(); // OJO, ver esto, el tamaño de ese vector debería ser el número de bloques

      P.resize(nBlocks, NULL);
      block_select.resize(nBlocks, NULL);
      block_rank.resize(nBlocks, NULL);

      B = sdsl::bit_vector(nBlocks, 0);

      uint64_t start = 0;
      uint64_t first_elem = 1;

      uint64_t last_elem, end, size_block, i = 0;
      for(i = 0; i < nBlocks - 1; i++){
        last_elem = partition[i];
        end = select_1(last_elem - 1);
        size_block = end - start + 1;

        block_bv.resize(size_block);
        sdsl::util::set_to_value(block_bv, 0);

        //last element is not include
        uint64_t amount_ones, elem, temp;
        for(elem = first_elem, amount_ones = 0; elem < last_elem; elem++, amount_ones++) {
          temp = select_1(elem);
          block_bv[temp - start] = 1;
        }
        elements_of_L.push_back(temp);
        elements_of_E.push_back(amount_ones + (i == 0 ? 0 : elements_of_E[i - 1]));


        uint64_t type_encoding_block = type_encoding_vigna(end - first_elem + 1, amount_ones);

        add_block(block_bv, type_encoding_block, i);
        start = end + 1;
        first_elem = partition[i];
      }
      // last block
      last_elem = partition[i];
      end = u;
      size_block = end - start;

      block_bv.resize(size_block);
      sdsl::util::set_to_value(block_bv, 0);

      uint64_t amount_ones, elem, temp;
      for(elem = first_elem, amount_ones = 0; elem <= last_elem; elem++, amount_ones++) {
        temp = select_1(elem);
        block_bv[temp - start] = 1;
      }
      elements_of_L.push_back(temp);
      elements_of_E.push_back(amount_ones + (nBlocks == 1 ? 0 : elements_of_E[i - 1]));
      uint64_t type_encoding_block = type_encoding_vigna(end - first_elem, amount_ones);

      add_block(block_bv, type_encoding_block, i);
      
      L = sux::bits::EliasFano<>(elements_of_L, elements_of_L[elements_of_L.size() - 1] + 1);
      E = sux::bits::EliasFano<>(elements_of_E, elements_of_E[elements_of_E.size() - 1] + 1);
      
      size_in_bits_auxiliar = size_in_bits_formula();
    } 

    pef_vector_opt_vigna(std::vector<uint64_t> &pb, uint64_t universe, double eps1, double eps2) {
      // setup fixed_cost 
      fixed_cost = _fixed_cost;

      sdsl::bit_vector block_bv;
      std::vector<uint64_t> elements_of_L;
      std::vector<uint64_t> elements_of_E;

      u = universe;
      n = pb.size();
  
      //std::cout << "ENTERING PEF VIGNA\n";
      //----- Optimal partition -----
      pair<std::vector<uint64_t>, uint64_t> ret = optimal_partition(pb, eps1, eps2);
      std::vector<uint64_t> partition = ret.first;
      uint64_t cost_opt = ret.second;
      //-----------------------------

      //std::cout << "BLOCKS SELECTED\n";
      nBlocks = partition.size(); // OJO, ver esto, el tamaño de ese vector debería ser el número de bloques
      cout << nBlocks << " " << cost_opt << " " << (double) cost_opt / n << "\n";
      //print(partition);

      P.resize(nBlocks, NULL);
      block_select.resize(nBlocks, NULL);
      block_rank.resize(nBlocks, NULL);

      B = sdsl::bit_vector(nBlocks, 0);

      uint64_t start = 0;
      uint64_t first_elem = 1;

      uint64_t last_elem, end, size_block, i = 0;
      for(i = 0; i < nBlocks - 1; i++){
        last_elem = partition[i];
        end = pb[last_elem - 1 - 1];
        size_block = end - start + 1;

        block_bv.resize(size_block);
        sdsl::util::set_to_value(block_bv, 0);

        //last element is not include
        uint64_t amount_ones, elem, temp;
        for(elem = first_elem, amount_ones = 0; elem < last_elem; elem++, amount_ones++) {
          temp = pb[elem - 1];
          block_bv[temp - start] = 1;
        }
        elements_of_L.push_back(temp);
        elements_of_E.push_back(amount_ones + (i == 0 ? 0 : elements_of_E[i - 1]));

        uint64_t type_encoding_block = type_encoding_vigna(end - start + 1, amount_ones);

        //std::cout << "ENTER FUNCTION ADD BLOCK\n";
        add_block(block_bv, type_encoding_block, i);
        start = end + 1;
        first_elem = partition[i];
      }
      //std::cout << "LAST BLOCK\n";
      // last block
      last_elem = partition[i];
      end = u;
      size_block = end - start;

      block_bv.resize(size_block);
      sdsl::util::set_to_value(block_bv, 0);

      uint64_t amount_ones, elem, temp;
      for(elem = first_elem, amount_ones = 0; elem <= last_elem; elem++, amount_ones++) {
        temp = pb[elem - 1];
        block_bv[temp - start] = 1;
      }
      elements_of_L.push_back(temp);
      elements_of_E.push_back(amount_ones + (nBlocks == 1 ? 0 : elements_of_E[i - 1]));
      
      uint64_t type_encoding_block = type_encoding_vigna(end - start, amount_ones);

      add_block(block_bv, type_encoding_block, i);

      L = sux::bits::EliasFano<>(elements_of_L, elements_of_L[elements_of_L.size() - 1] + 1);

      E = sux::bits::EliasFano<>(elements_of_E, elements_of_E[elements_of_E.size() - 1] + 1);

      size_in_bits_auxiliar = size_in_bits_formula();
      //std::cout << "FINAL PEF OPT\n";
    } 

    uint64_t select(uint64_t i){
      assert(i <= n);
      
      uint64_t blk = E.rank(i);
      if (B[blk]) {
        if (blk == 0) {
          return  (*(sux::bits::EliasFano<> *)P[blk]).select(i - 1);
        } else {
          uint64_t i_block = i - E.select(blk - 1);
          return L.select(blk - 1) + 1 + (*(sux::bits::EliasFano<> *)P[blk]).select(i_block - 1); 
        }
      } else {
        if (P[blk]) {
          if (blk == 0) {
              
            return (*(select_support_c *)block_select[blk])(i);
          } else {
            uint64_t i_block = i - E.select(blk - 1);
            return L.select(blk - 1) + 1 + (*(select_support_c *)block_select[blk])(i_block);
          }
        } else {
          // block is a run of 1s
          if (blk == 0) {
            return (i - 1);
          } else {            
            uint64_t i_block = i - E.select(blk - 1);
            return L.select(blk - 1) + i_block;
          }
        }
      }

      return 0;
    }

    uint64_t rank(uint64_t i){
      if(i >= u) return n_ones();

      uint64_t blk = L.rank(i);
      if(blk == nBlocks) return n_ones();

      uint64_t rank_val = 0;

      if(blk > 0)
        rank_val += E.select(blk - 1);

      if (B[blk]) {
        if(blk == 0){
          rank_val += (*(sux::bits::EliasFano<> *)P[blk]).rank(i);
        } else {
          //obtain position in the block
          uint64_t i_block = i - 1 - L.select(blk - 1);
          rank_val += (*(sux::bits::EliasFano<> *)P[blk]).rank(i_block);
        }
      } else {
        if (P[blk]) {
          if(blk == 0){
            rank_val += (*(rank_support_c *)block_rank[blk])(i);
          } else{
            
            //obtain position in the block
            uint64_t i_block = i - 1 - L.select(blk - 1);
            rank_val += (*(rank_support_c *)block_rank[blk])(i_block);
          }
        } else {
          // block is a run of 1s
          if(blk == 0){
            rank_val += i;
          } else{
            rank_val += i - 1 - L.select(blk - 1);
          }
        }
      }
      return rank_val;
    }

    uint64_t operator[](uint64_t i){
      //if i >= u doesn't make sense

      uint64_t amount_ones = rank(i + 1);
      if(amount_ones == 0) return 0;

      uint64_t pos_last_one = select(amount_ones);

      if(pos_last_one < i) return 0;
      else return 1; 
    }
  
};

#endif
