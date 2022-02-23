#ifndef PEF_UNIF
#define PEF_UNIF

#include <vector>
#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <cinttypes>
#include <assert.h> 

using namespace std;
using namespace sdsl;


template<uint64_t b=256>
class pef_vector_unif {
  sdsl::sd_vector<> L;
  sdsl::select_support_sd<1> select_L;
  sdsl::rank_support_sd<1> rank_L;


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
    
public:

  pef_vector_unif() {;}

  ~pef_vector_unif() {;}

  uint64_t size() {
    return u;
  }

  uint64_t n_ones() {
    return n;
  }

  uint64_t n_blocks() {
    return nBlocks;
  }

  uint64_t size_in_bytes() {
    uint64_t size = sdsl::size_in_bytes(L) + sdsl::size_in_bytes(select_L) 
                  + sdsl::size_in_bytes(rank_L) 
                  + sdsl::size_in_bytes(B) + 3 * sizeof(uint64_t) 
                  + nBlocks*sizeof(void *);
        
    for (uint64_t i=0; i < nBlocks; ++i) {
      if (B[i]) {
        size += sdsl::size_in_bytes(*(sd_vector<> *)P[i]) 
              + sdsl::size_in_bytes(*(select_support_sd<1> *)block_select[i]);
              + sdsl::size_in_bytes(*(rank_support_sd<1> *)block_rank[i]);
      } else {
        if (P[i]) {
          size += sdsl::size_in_bytes(*(bit_vector *)P[i]) 
                + sdsl::size_in_bytes(*(select_support_scan<1> *)block_select[i])
                + sdsl::size_in_bytes(*(rank_support_scan<1> *)block_rank[i]);
        }
      }
    }
    return size;
  }

  pef_vector_unif(bit_vector &bv) {
    
    select_support_mcl<1> select_1(&bv);
    rank_support_v<1> rank_1(&bv);
    
    bit_vector block_bv;
    std::vector<uint64_t> elements_of_L;

    u = bv.size();
    n = rank_1(bv.size());
    nBlocks = (n + b - 1) / b;
          
          
    P.resize(nBlocks, NULL);
    block_select.resize(nBlocks, NULL); 
    block_rank.resize(nBlocks, NULL); 

    B = bit_vector(nBlocks, 0);

    uint64_t next_elem, nElems, start = 0, univ_size_block, i, temp;

    //cout << "Number of 1s = " << n << " Number of blocks = " << nBlocks << endl;
    for (i = 0; i < nBlocks - 1; i++) {
      //cout << "last 1 in block " << select_1((i+1)*b) << endl; 
      univ_size_block = select_1((i+1)*b) - start + 1;

      //cout << "Block " << i << " universe size " << univ_size_block << endl;
      block_bv.resize(univ_size_block);
      sdsl::util::set_to_value(block_bv, 0);

      for (next_elem = i*b+1, nElems = 0; nElems < b; ++nElems, ++next_elem) {
        temp = select_1(next_elem);
        block_bv[temp-start] = 1;
      }
      elements_of_L.push_back(temp);
      // Now, chooses the encoding of the block
      if (univ_size_block == b) {
        B[i] = 0;
        P[i] = NULL; 
        block_select[i] = NULL;
        block_rank[i] = NULL;
      } else if (b > univ_size_block/4) {
        B[i] = 0;
        P[i] = new bit_vector(block_bv);
        block_select[i] = new select_support_scan<1>((bit_vector *)P[i]); 
        block_rank[i] = new rank_support_scan<1>((bit_vector *)P[i]); 
      } else {
        B[i] = 1;
        P[i] = new sd_vector<>(block_bv);
        block_select[i] = new select_support_sd<1>((sd_vector<> *)P[i]);
        block_rank[i] = new rank_support_sd<1>((sd_vector<> *)P[i]);
      }

      //cout << "next elem " << next_elem << endl;
      start = select_1(next_elem-1)+1;
      //cout << "start " << start << endl;
    }

    // last block
    univ_size_block = u - start;  // revisar esta formula
    //cout << "Block " << i << " universe size " << univ_size_block << endl;

    block_bv.resize(univ_size_block);
    sdsl::util::set_to_value(block_bv, 0);

    for(next_elem = i*b+1; next_elem <= n; next_elem++) {
      temp = select_1(next_elem);
      block_bv[temp-start] = 1;
    }

    elements_of_L.push_back(temp);

    if (univ_size_block == b) {
      B[i] = 0;
      P[i] = NULL; 
      block_select[i] = NULL;
      block_rank[i] = NULL;
    } else if (b > univ_size_block / 4) {
      B[i] = 0;
      P[i] = new bit_vector(block_bv);
      block_select[i] = new select_support_scan<1>((bit_vector *)P[i]); 
      block_rank[i] = new rank_support_scan<1>((bit_vector *)P[i]); 
    } else {
      B[i] = 1;
      P[i] = new sd_vector<>(block_bv);
      block_select[i] = new select_support_sd<1>((sd_vector<> *)P[i]);
      block_rank[i] = new rank_support_sd<1>((sd_vector<> *)P[i]);
    }


    L = sd_vector<>(elements_of_L.begin(), elements_of_L.end());
    util::init_support(select_L, &L);
    util::init_support(rank_L, &L);
  }


  uint64_t select(uint64_t i) {
    assert(i <= n);
    uint64_t blk = (i - 1) / b;
    if (B[blk]) {
      if (blk == 0) {
        return (*(select_support_sd<1> *)block_select[blk])(i);
      } else {
        //cout << "PEF" << endl; fflush(stdout);
        return select_L(blk) + 1 + (*(select_support_sd<1> *)block_select[blk])((i - 1) % b + 1); 
      }
    } else {
      if (P[blk]) {
        if (blk == 0) {
          return (*(select_support_scan<1> *)block_select[blk])(i);
        } else {
          //cout << "Plain" << endl; fflush(stdout);
          return select_L(blk) + 1 + (*(select_support_scan<1> *)block_select[blk])((i - 1) % b + 1);
        }
      } else {
        // block is a run of 1s
        if (blk == 0) {
          return (i - 1);
        } else {
          return select_L(blk) + (i - 1) % b + 1;
        }
      }
    }
  }

  uint64_t rank(uint64_t i) {
    if(i >= u) return n_ones();

    uint64_t blk = rank_L(i); //obtain the block of the i bit

    if(blk == nBlocks) return n_ones();
    
    uint64_t rank_val = 0;

    if(blk > 0)
      rank_val += (b * blk); //for each block there is b 1's
                             //rank value until the block is blk * b

    if (B[blk]) {
      if(blk == 0){
        rank_val += (*(rank_support_sd<1> *)block_rank[blk])(i);
      } else {
        //obtain position in the block
        uint64_t i_block = i - 1 - select_L(blk);
        rank_val += (*(rank_support_sd<1> *)block_rank[blk])(i_block);
      }
    } else {
      if (P[blk]) {
        if(blk == 0){
          rank_val += (*(rank_support_scan<1> *)block_rank[blk])(i);
        } else{
          //obtain position in the block
          uint64_t i_block = i - 1 - select_L(blk);
          rank_val += (*(rank_support_scan<1> *)block_rank[blk])(i_block);
        }
      } else {
        // block is a run of 1s
        if(blk == 0){
          rank_val += i;
        } else{
          rank_val += i - 1 - select_L(blk);
        }
      }
    }
    

    return rank_val;
  }  
};

#endif
