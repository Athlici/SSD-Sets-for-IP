#include <vector>
#include <future>
#include <fstream>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "BS_thread_pool.hpp" //https://github.com/bshoshany/thread-pool
#include "boost/dynamic_bitset.hpp"

using std::vector;
using std::future;
using Eigen::VectorXi;
using boost::dynamic_bitset;

//Take the absolute value of an integer
int abs(int a){
  return (a>=0) ? a : -a;
}

//Generates all bounded vectors of dimension m.
//The bound c is on Delta for a=0, on |A|1 for a=1.
vector<VectorXi> vecs(int a, int m, int c){
  vector<VectorXi> res = {};
  VectorXi tmp = VectorXi::Zero(m);

  auto vecGen = [a,m,res,tmp](int k, int r, bool z, auto& vGref) mutable -> vector<VectorXi> {
    if(k == m)
      res.push_back(tmp);
    else
      for(tmp(k)=r*(z-1); tmp(k)<=r; tmp(k)++)
        vGref(k+1, r-a*abs(tmp(k)), z && tmp(k)==0, vGref);
    return res;
  };

  return vecGen(0, c, true, vecGen);
}

vector<vector<VectorXi>> setFind(uint8_t s, uint16_t c, vector<VectorXi> v, uint16_t s0, uint16_t s1){
  uint8_t m = 2;

  vector<vector<VectorXi>> res = {};

  uint16_t S[s+1] = {s0,s1}; S[2] = s1;
  uint16_t max[s+1][m];
  uint16_t min[s+1][m];
  for(int i=0; i<m; i++){
    max[1][i] = min[1][i] = 0;
    if(v[s0](i) >= 0)
      max[1][i] += v[s0](i);
    else
      min[1][i] -= v[s0](i);
    if(v[s1](i) >= 0)
      max[1][i] += v[s1](i);
    else
      min[1][i] -= v[s1](i);
  }

  //VectorXi cord = VectorXi::Ones(m), off = VectorXi::Ones(m);
  int d = 2*s*c+1;
  //cord(0) = 0; cord(1) = 1;
  //for(int i=2; i<m; i++)
  //  cord(i) = cord(i-1)*d;
  int off = s*c;
  //d *= cord(m-1);

  dynamic_bitset<> p[s][d];
  for(int i=1; i<s; i++)
    for(int j=0; j<d; j++){
      dynamic_bitset<> pij(s*c);
      p[i][j] = pij;
    }
  p[1][off][0] = 1;
  p[1][off+v[s0](1)][v[s0](0)] = 1;
  p[1][off+v[s1](1)][v[s1](0)] = 1;
  p[1][off+v[s0](1)+v[s1](1)][v[s0](0)+v[s1](0)] = 1;

  for(uint8_t k=1; k>0; k--){
    uint16_t Sk = S[k];
    for(uint16_t Sk1 = S[k+1]+1; Sk1 < v.size(); Sk1++){
      if(!p[k][off+v[Sk1](1)][v[Sk1](0)]){
        bool broke = false;
        for(int16_t i = off-min[k][1]; i<=off+max[k][1]; i++){
          int ind = i+v[Sk1](1);
          p[k+1][ind] = p[k][i];
          p[k+1][ind] <<= v[Sk1](0);
          p[k+1][ind] &= p[k][ind];
          if(p[k+1][ind].any()){
            broke = true;
            break;
          }
        }
      if(!broke){
        S[k+1] = Sk1;
        if(k+2<s){
          if(v[Sk1](1)>0){
            min[k+1][1] = min[k][1];
            max[k+1][1] = max[k][1]+v[Sk1](1);
          } else {
            min[k+1][1] = min[k][1]-v[Sk1](1);
            max[k+1][1] = max[k][1];
          }

          for(int16_t i = off-min[k][1]; i<=off+max[k][1]; i++)
            p[k+1][i] = p[k][i];
          for(int16_t i = off-min[k][1]; i<=off+max[k][1]; i++){
            int ind = i+v[Sk1](1);
            p[k+1][ind] |= (p[k][i] << v[Sk1](0));
          }

          S[k+2] = Sk1;
          k += 2;
        } else {
          std::cout << "{";
          for(uint8_t i=0; i<s; i++){
            std::cout << "{";
            for(uint8_t j=0; j<m; j++)
              std::cout << v[S[i]](j) << ",";
            std::cout << "\b},";
          }
          std::cout << "\b}\n";
          vector<VectorXi> Sv(s);
          for(uint8_t i=0; i<s; i++)
            Sv.push_back(v[S[i]]);
          res.push_back(Sv);
          k += 1;
        }
        break;
      }}
    }
  }

  return res;
}

int main(int argc, char** argv){
  //Take 4 arguments from the command line
  uint8_t a = atoi(argv[1]);    //0 -> Delta, 1 -> |A|1 norm for bounds
  uint8_t m = 2;    //dimension
  uint8_t s = atoi(argv[2]);    //desired set size
  uint16_t l = atoi(argv[3]), u = atoi(argv[4]); //lower and upper bound

  BS::thread_pool pool;
  vector<future<vector<vector<VectorXi>>>> futes = {};

  uint16_t ind[u+1]; ind[l-1]=0;
  for(uint16_t k=l; k<=u; k++){
    vector<VectorXi> v = vecs(a,2,k);
    uint16_t vs = v.size(), vts = v.size();
    for(uint16_t i=1; i<vts; i++)
      for(uint16_t j=i+1; j<vs; j++)
        //TODO: remove symmetry duplicates for speedup
        futes.push_back(pool.submit(setFind, s, k, v, i, j));
    ind[k] = futes.size();
  }

  for(uint16_t k=l; k<=u; k++){
    vector<vector<VectorXi>> res = {};
    for(uint16_t j=ind[k-1]; j<ind[k]; j++){
      vector<vector<VectorXi>> resA = futes[j].get();
      res.insert(res.end(), resA.begin(), resA.end());
    }

    std::ofstream out;
    out.open("out"+std::to_string(k)+".txt");
    out << 2 << " " << +s << "\n";
    for(auto S : res){
      for(auto v : S)
        for(int i=0; i<v.size(); i++)
          out << v(i) << " ";
      out << "\n";
    }
    out.close();
    std::cout << "Finished searching: " << k << "\n";
  }
}
