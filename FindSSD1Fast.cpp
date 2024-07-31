#include <vector>
#include <future>
#include <fstream>
#include <iostream>
#include "BS_thread_pool.hpp" //https://github.com/bshoshany/thread-pool

using std::vector;
using std::future;

//use known values to do minor pruning, subset of SSD is SSD.
const uint16_t smin[] = {1,2,4,7,13,24,44,84,161,309};

//This function is effectively a tuned implementation of
// return (pk & (pk << Sk)) > 0;
//for the bitset pk of length mk shifted by Sk.
bool shiftAndAny(uint16_t mk, uint16_t Sk, uint64_t* pk){
  uint16_t maxk = mk/64;        //largest index
  uint16_t offset = Sk/64;      //inter bit shift
  uint16_t shift1 = Sk%64;      //intra bit shift
  uint16_t shift2 = 64-shift1;
  uint64_t a, b = 0, c;
  // x << 64 does not change x, instead of zeroing it out.
  //Therefore we need a case distinction for Sk divisible by 64
  if(shift1>0){
    //When the offsets are small we can go faster with few variables
    switch(offset){
      case 0:
        for(uint16_t i = 0; i <= maxk; i++){
          a = pk[i];
          b += a << shift1;
          if(a & b)
            return true;
          b = a >> shift2;
        }
        break;
      case 1:
        c = pk[0];
        for(uint16_t i = 1; i <= maxk; i++){
          a = pk[i];
          b += c << shift1;
          if(a & b)
            return true;
          b = c >> shift2;
          c = a;
        }
        break;
      default:
        for(uint16_t i = offset; i <= maxk; i++){
          a = pk[i];
          c = pk[i-offset];
          b += c << shift1;
          if(a & b)
            return true;
          b = c >> shift2;
        }
    }
  } else {
    for(uint16_t i = offset; i <= maxk; i++)
      if(pk[i] & pk[i-offset])
        return true;
  }
  return false;
}

//This function is effectively a tuned implementation of
// pk1 = pk | (pk << Sk);
//for bitsets pk and pk1 of length mk shifted by Sk
void shiftOr(uint16_t mk, uint16_t Sk, uint64_t* pk, uint64_t* pk1){
  uint16_t maxk = mk/64;
  uint16_t maxk1 = (mk+Sk)/64;
  uint16_t offset = Sk/64;
  uint16_t shift1 = Sk%64;
  uint16_t shift2 = 64-shift1;
  uint64_t a, b = 0, c;
  //Again the divisible by 64 distinction
  if(shift1>0){
    //Specialisation for Sk<64,128 did not increase speed, hence omitted
    b = 0;
    for(uint16_t i = 0; i < offset; i++)
      pk1[i] = pk[i];
    for(uint16_t i = offset; i <= maxk; i++){
      c = pk[i-offset];
      b += c << shift1;
      pk1[i] = pk[i] | b;
      if(shift2==64)
        b = 0;
      else
      b = c >> shift2;
    }
    for(uint16_t i = maxk+1; i <= maxk1; i++){
      c = pk[i-offset];
      b += c << shift1;
      pk1[i] = b;
      b = c >> shift2;
    }
  } else {
    for(uint16_t i = 0; i < offset; i++)
      pk1[i] = pk[i];
    for(uint16_t i = offset; i <= maxk; i++)
      pk1[i] = pk[i] | pk[i-offset];
    for(uint16_t i = maxk+1; i <= maxk1; i++)
      pk1[i] = pk[i-offset];
  }
}

//Find ssd sets of size s with largest entries s0 and s1
vector<vector<uint16_t>> setFind(uint8_t s, uint16_t s0, uint16_t s1){
  vector<vector<uint16_t>> res = {};

  uint16_t S[s] = {s0,s1}; S[2] = s1;
  uint16_t max[s]; max[1] = s0+s1;
  uint64_t*  p[s];
  for(int i=1; i<s; i++)
    p[i] = (uint64_t*) calloc(1+s0*s/64,8);
  p[1][0]          |= 1ULL;
  p[1][s0/64]      |= 1ULL << (s0%64);
  p[1][s1/64]      |= 1ULL << (s1%64);
  p[1][(s0+s1)/64] |= 1ULL << ((s0+s1)%64);

  //Perform depth first search for SSD sets
  for(uint8_t k=1; k>0; k--){
    uint16_t Sk = S[k], Smin = smin[s-k-2];
    for(uint16_t Sk1 = S[k+1]-1; Sk1 >= Smin; Sk1--){
      if(!((p[k][Sk1/64]>>(Sk1%64))&1) && //Only take numbers not yet reached
         !shiftAndAny(max[k],Sk1,p[k])){  //which don't create duplicates
        S[k+1] = Sk1;
        if(k+2<s){  //Add another number to the current set
          max[k+1] = max[k]+Sk1;
          p[k+1][1+max[k+1]/64] = 0;
          shiftOr(max[k],Sk1,p[k],p[k+1]);
          S[k+2] = Sk1;
          k += 2;
        } else {    //Or print and store the result
          std::cout << "{";
          for(uint8_t i=0; i<s; i++)
            std::cout << +S[i] << ",";
          std::cout << "\b}\n";
          vector<uint16_t> Sv(S,S+s);
          res.push_back(Sv);
          k += 1;
        }
        break;
      }
    }
  }

  return res;
}

int main(int argc, char** argv){

  //Take 3 arguments from the command line, the set size, lower and upper bound of values
  uint8_t s = atoi(argv[1]);
  uint16_t l = atoi(argv[2]), u = atoi(argv[3]);

  BS::thread_pool pool;     //Automatically uses the available number of cores

  vector<future<vector<vector<uint16_t>>>> futes = {};
  for(uint16_t i=l; i<=u; i++)
    for(uint16_t j=i-1; j>=smin[s-2]; j--)
      futes.push_back(pool.submit(setFind, s, i, j));

  int k=0;
  for(uint16_t i=l; i<=u; i++){
    vector<vector<uint16_t>> res = {};
    for(uint16_t j=i-1; j>=smin[s-2]; j--){
      vector<vector<uint16_t>> resA = futes[k].get();
      res.insert(res.end(), resA.begin(), resA.end());
      k += 1;
    }

    //Write all found sets to a file, output is unreliable with multithreading
    //(Could be fixed, BS offers some functionality, but not worth the effort)
    std::ofstream out;
    out.open("out"+std::to_string(i)+".txt");
    out << 1 << " " << +s << "\n";
    for(auto S : res){
      for(auto v : S)
        out << +v << " ";
      out << "\n";
    }
    out.close();
    std::cout << "Finished searching: " << i << "\n";
  }
}
