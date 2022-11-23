#include "pef_vector_unif.hpp"
#include "pef_vector_opt.hpp"
#include "pef_vector_opt_vigna.hpp"
#include "pef_vector_opt_il.hpp"

//#define BITS 20000000
//#define N 196433254820

int main(int argc, char** argv){
  vector< uint64_t > pb;
  // argv will contain all the test files
  // need to be with the format `universe filepath` 
  for(uint64_t i = 1; i < argc; i += 2){
    uint64_t universe = strtoull(argv[i], NULL, 0);
    char* filepath = argv[i + 1];
    cout << universe << " " << filepath << "\n";
    
    ifstream rf(filepath, std::ios::binary);
    if(!rf) {
      std::cout << "Cannot open file! " << filepath << "\n";
      continue;
    }

    uint64_t pos_one;
    uint64_t readed = 0;
    
    while(!rf.eof()){
      rf.read((char*) &pos_one, sizeof(uint64_t));
      pb.push_back(pos_one);
      //readed++;
      //if(readed >= BITS) break;
    }
    rf.close();
  

    //for(uint64_t tests = 0; tests < 200; tests++){
    //  sdsl::bit_vector bv(BITS);
    //  sdsl::util::set_random_bits(bv);
    //  for(uint64_t i = 0; i < BITS; i++){
    //    if(bv[i]) pb.push_back(i);
    //  }

    //universe = pb[pb.size() - 1] + 1; 
/*
    sdsl::bit_vector bv(universe, 0);
    for(auto bit : pb) bv[bit] = 1;
    
    sdsl::rank_support_v<1> rank_bv(&bv);
    sdsl::select_support_mcl<1> select_bv(&bv);
    cout << "SIZE BIT VECTOR " 
         << (double) (sdsl::size_in_bytes(bv) + sdsl::size_in_bytes(rank_bv) + sdsl::size_in_bytes(select_bv)) * 8 / pb.size()
         << "\n";

*/
    sdsl::sd_vector<> sd(pb.begin(), pb.end());
    sdsl::rank_support_sd<1> rank_sd(&sd);
    sdsl::select_support_sd<1> select_sd(&sd);
    cout << "SIZE SD VECTOR " 
         << (double) (sdsl::size_in_bytes(sd) + sdsl::size_in_bytes(rank_sd) + sdsl::size_in_bytes(select_sd)) * 8 / pb.size()
         << "\n";
    //cout << "PLAIN_BIT_VECTOR FORMULA " << (double) bitsize_plain_bitvector(universe, BITS) / pb.size() << "\n";
    //cout << "SD_VECTOR FORMULA " << (double) bitsize_elias_fano(universe, BITS) / pb.size() << "\n";
    //pef_vector_unif<4096, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> pef_unif(pb, universe);
    //cout << "SIZE PEF_UNIF " << ((double)pef_unif.size_in_bits_formula()) / ((double) pb.size()) << "\n"; 

    pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 64> pef_opt(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT " << pef_opt.size_in_bits_formula() << " " << ((double)pef_opt.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 128> pef_opt2(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT " << pef_opt2.size_in_bits_formula() << " " << ((double)pef_opt2.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 256> pef_opt3(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT " << pef_opt3.size_in_bits_formula() << " " << ((double)pef_opt3.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 512> pef_opt4(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT " << pef_opt4.size_in_bits_formula() << " " << ((double)pef_opt4.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 1024> pef_opt5(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT " << pef_opt5.size_in_bits_formula() << " " << ((double)pef_opt5.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
 
    pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 2048> pef_opt6(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT " << pef_opt6.size_in_bits_formula() << " " << ((double)pef_opt6.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
  
    pef_vector_opt_vigna<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 64> pef_opt_vigna(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT VIGNA " << pef_opt_vigna.size_in_bits_formula() << " " << ((double)pef_opt_vigna.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt_vigna<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 128> pef_opt2_vigna(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT VIGNA " << pef_opt2_vigna.size_in_bits_formula() << " " << ((double)pef_opt2_vigna.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt_vigna<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 256> pef_opt3_vigna(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT VIGNA " << pef_opt3_vigna.size_in_bits_formula() << " " << ((double)pef_opt3_vigna.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt_vigna<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 512> pef_opt4_vigna(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT VIGNA " << pef_opt4_vigna.size_in_bits_formula() << " " << ((double)pef_opt4_vigna.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
  
    pef_vector_opt_vigna<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 1024> pef_opt5_vigna(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT VIGNA " << pef_opt5_vigna.size_in_bits_formula() << " " << ((double)pef_opt5_vigna.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt_vigna<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 2048> pef_opt6_vigna(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT VIGNA " << pef_opt6_vigna.size_in_bits_formula() << " " << ((double)pef_opt6_vigna.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
/*
    pef_vector_opt_il<sdsl::rank_support_il<1, 128>, sdsl::select_support_il<1, 128>, 1, 64, 128> pef_opt_il_128_1(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 128 " << ((double)pef_opt_il_128_1.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 128>, sdsl::select_support_il<1, 128>, 1, 128, 128> pef_opt_il_128_2(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 128 " << ((double)pef_opt_il_128_2.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt_il<sdsl::rank_support_il<1, 128>, sdsl::select_support_il<1, 128>, 1, 256, 128> pef_opt_il_128_3(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 128 " << ((double)pef_opt_il_128_3.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt_il<sdsl::rank_support_il<1, 128>, sdsl::select_support_il<1, 128>, 1, 512, 128> pef_opt_il_128_4(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 128 " << ((double)pef_opt_il_128_4.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt_il<sdsl::rank_support_il<1, 128>, sdsl::select_support_il<1, 128>, 1, 1024, 128> pef_opt_il_128_5(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 128 " << ((double)pef_opt_il_128_5.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
    pef_vector_opt_il<sdsl::rank_support_il<1, 128>, sdsl::select_support_il<1, 128>, 1, 2048, 128> pef_opt_il_128_6(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 128 " << ((double)pef_opt_il_128_6.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 256>, sdsl::select_support_il<1, 256>, 1, 64, 256> pef_opt_il_256_1(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 256 " << ((double)pef_opt_il_256_1.size_in_bits_formula()) / ((double) pb.size()) << "\n";

    pef_vector_opt_il<sdsl::rank_support_il<1, 256>, sdsl::select_support_il<1, 256>, 1, 128, 256> pef_opt_il_256_2(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 256 " << ((double)pef_opt_il_256_2.size_in_bits_formula()) / ((double) pb.size()) << "\n";
    
    pef_vector_opt_il<sdsl::rank_support_il<1, 256>, sdsl::select_support_il<1, 256>, 1, 256, 256> pef_opt_il_256_3(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 256 " << ((double)pef_opt_il_256_3.size_in_bits_formula()) / ((double) pb.size()) << "\n";
    
    pef_vector_opt_il<sdsl::rank_support_il<1, 256>, sdsl::select_support_il<1, 256>, 1, 512, 256> pef_opt_il_256_4(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 256 " << ((double)pef_opt_il_256_4.size_in_bits_formula()) / ((double) pb.size()) << "\n";
    
    pef_vector_opt_il<sdsl::rank_support_il<1, 256>, sdsl::select_support_il<1, 256>, 1, 1024, 256> pef_opt_il_256_5(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 256 " << ((double)pef_opt_il_256_5.size_in_bits_formula()) / ((double) pb.size()) << "\n";
    
    pef_vector_opt_il<sdsl::rank_support_il<1, 256>, sdsl::select_support_il<1, 256>, 1, 2048, 256> pef_opt_il_256_6(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 256 " << ((double)pef_opt_il_256_6.size_in_bits_formula()) / ((double) pb.size()) << "\n";

    pef_vector_opt_il<sdsl::rank_support_il<1, 512>, sdsl::select_support_il<1, 512>, 1, 64, 512> pef_opt_il_512_1(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 512 " << ((double)pef_opt_il_512_1.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 512>, sdsl::select_support_il<1, 512>, 1, 128, 512> pef_opt_il_512_2(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 512 " << ((double)pef_opt_il_512_2.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 512>, sdsl::select_support_il<1, 512>, 1, 256, 512> pef_opt_il_512_3(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 512 " << ((double)pef_opt_il_512_3.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 512>, sdsl::select_support_il<1, 512>, 1, 512, 512> pef_opt_il_512_4(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 512 " << ((double)pef_opt_il_512_4.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 512>, sdsl::select_support_il<1, 512>, 1, 1024, 512> pef_opt_il_512_5(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 512 " << ((double)pef_opt_il_512_5.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 512>, sdsl::select_support_il<1, 512>, 1, 2048, 512> pef_opt_il_512_6(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 512 " << ((double)pef_opt_il_512_6.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 1024>, sdsl::select_support_il<1, 1024>, 1, 64, 1024> pef_opt_il_1024_1(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 1024 " << ((double)pef_opt_il_1024_1.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 1024>, sdsl::select_support_il<1, 1024>, 1, 128, 1024> pef_opt_il_1024_2(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 1024 " << ((double)pef_opt_il_1024_2.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 1024>, sdsl::select_support_il<1, 1024>, 1, 256, 1024> pef_opt_il_1024_3(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 1024 " << ((double)pef_opt_il_1024_3.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 1024>, sdsl::select_support_il<1, 1024>, 1, 512, 1024> pef_opt_il_1024_4(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 1024 " << ((double)pef_opt_il_1024_4.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 1024>, sdsl::select_support_il<1, 1024>, 1, 1024, 1024> pef_opt_il_1024_5(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 1024 " << ((double)pef_opt_il_1024_5.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 1024>, sdsl::select_support_il<1, 1024>, 1, 2048, 1024> pef_opt_il_1024_6(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 1024 " << ((double)pef_opt_il_1024_6.size_in_bits_formula()) / ((double) pb.size()) << "\n";  


 //     pef_vector_opt_vigna<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>, 2, 192> pef_opt_vigna(pb, universe, 0.03, 0.3);
   //   std::cout << "SIZE PEF_OPT_VIGNA " << ((double)pef_opt_vigna.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
  
    pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 1, 1, 192> pef_opt_XD(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT NO SR " << ((double)pef_opt_XD.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
      pef_vector_opt_vigna<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>, 1, 192> pef_opt_vigna_XD(pb, universe, 0.03, 0.3);
      std::cout << "SIZE PEF_OPT_VIGNA NO SR " << ((double)pef_opt_vigna_XD.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 2, 2, 96> pef_opt_XD2(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT NO SR " << ((double)pef_opt_XD2.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
      pef_vector_opt_vigna<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>, 1, 96> pef_opt_vigna_XD2(pb, universe, 0.03, 0.3);
      std::cout << "SIZE PEF_OPT_VIGNA NO SR " << ((double)pef_opt_vigna_XD2.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
pef_vector_opt<sdsl::rank_support_v<1>, sdsl::select_support_mcl<1>, 2, 2, 64> pef_opt_XD3(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT NO SR " << ((double)pef_opt_XD3.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    
      pef_vector_opt_vigna<sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>, 2, 192> pef_opt_vigna_XD3(pb, universe, 0.03, 0.3);
      std::cout << "SIZE PEF_OPT_VIGNA NO SR " << ((double)pef_opt_vigna_XD3.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    */
      //pb.clear();
      //sdsl::util::clear(bv);
    //}
    
    /*
    pef_vector_opt_il<sdsl::rank_support_il<1, 128>, sdsl::select_support_il<1, 128>, 1, 1024, 128> pef_opt_il_128(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 128 " << ((double)pef_opt_il_128.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 256>, sdsl::select_support_il<1, 256>, 1, 1024, 256> pef_opt_il_256(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 256 " << ((double)pef_opt_il_256.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 512>, sdsl::select_support_il<1, 512>, 1, 1024, 512> pef_opt_il_512(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 512 " << ((double)pef_opt_il_512.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 1024>, sdsl::select_support_il<1, 1024>, 1, 1024, 1024> pef_opt_il_1024(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 1024 " << ((double)pef_opt_il_1024.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 2048>, sdsl::select_support_il<1, 2048>, 1, 1024, 2048> pef_opt_il_2048(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 2048 " << ((double)pef_opt_il_2048.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    pef_vector_opt_il<sdsl::rank_support_il<1, 4096>, sdsl::select_support_il<1, 4096>, 1, 1024, 4096> pef_opt_il_4096(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 4096 " << ((double)pef_opt_il_4096.size_in_bits_formula()) / ((double) pb.size()) << "\n";  
    pef_vector_opt_il<sdsl::rank_support_il<1, 4096>, sdsl::select_support_il<1, 4096>, 1, 2048, 4096> pef_opt_il_4096_fix(pb, universe, 0.03, 0.3);
    cout << "SIZE PEF_OPT_VIGNA_IL 4096 " << ((double)pef_opt_il_4096_fix.size_in_bits_formula()) / ((double) pb.size()) << "\n";  

    */


  }
  return 0;
}
