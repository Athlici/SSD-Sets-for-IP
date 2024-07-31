#include <vector>
#include <future>
#include <fstream>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "BS_thread_pool.hpp" //https://github.com/bshoshany/thread-pool

using std::vector;
using std::future;
using Eigen::MatrixXi;
using Eigen::VectorXi;

//Integer absolute value function
int abs(int a){
  return (a>=0) ? a : -a;
}

//Generates all norm bounded vectors of dimension m.
//The bound c is on Delta for a=0, on Amax for a=1.
vector<VectorXi> vecs(int a, int m, int c){
  vector<VectorXi> res = {};
  VectorXi tmp = VectorXi::Zero(m);

  //Recursive lambda function,
  //k is recursion depth, r is remaining norm, z is to only have v and not -v
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

//Finds subset sum distinct sets of size s, dimension m and bounded by c.
//ToDo: Eliminate redundancy of m! in RowPermutations?
vector<vector<VectorXi>> MaxSupp(int a, int m, int s, int g, int c){
  vector<VectorXi> v = vecs(a,m,c);
  sort(v.begin(), v.end(),
    [](const VectorXi & a, const VectorXi & b) -> bool
    { return a.lpNorm<1>() < b.lpNorm<1>(); });
  for(int j=0; j<m; j++){
    for(int i=1; i<v.size(); i++)
      std::cout << v[i](j) << " ";
    std::cout << "\n";
  }

  VectorXi cord = VectorXi::Ones(m), off = VectorXi::Ones(m);
  int d = 2*s*c+1;
  for(int i=1; i<m; i++)
    cord(i) = cord(i-1)*d;
  off *= s*c;

  auto subsetTest = [s,m,v,cord,off](int k, vector<VectorXi> sub, int l,
           vector<VectorXi> values, vector<bool> points){
  auto subsetTestSub = [s,m,v,cord,off](int k, vector<VectorXi> sub, int l,
           vector<VectorXi> values, vector<bool> points, auto& sTref){

    vector<vector<VectorXi>> res = {};
    if(k==s){
        std::cout << "{";
        for(auto vs : sub){
          std::cout << "{";
          for(int i=0; i<vs.size(); i++)
            std::cout << vs(i) << ",";
          std::cout << "\b},";
          }
        std::cout << "\b}\n";

      res.push_back(sub);
      return res;
    }

    for(int j=0; j<l; j++){
      VectorXi vj = v[j];
      if(!points[cord.dot(vj+off)]){
        bool broke = false;
        vector<    bool> pointNew = points;
        vector<VectorXi> valueNew = values;
        sub.push_back(vj);

        for(VectorXi a : values){
          VectorXi tmp = vj+a;
          int pos = cord.dot(tmp+off);
          if(pointNew[pos]){
            broke = true;
            break;
          } else {
            pointNew[pos] = true;
            valueNew.push_back(tmp);
          }
        }

        if(!broke){
          vector<vector<VectorXi>> resA = sTref(k+1, sub, j+1, valueNew, pointNew, sTref);
          res.insert(res.end(), resA.begin(), resA.end());
        }
        sub.pop_back();
      }
    }
    return res;
    };
    return subsetTestSub(k,sub,l,values,points,subsetTestSub);
  };

  BS::thread_pool pool;

  vector<future<vector<vector<VectorXi>>>> futes = {};
  VectorXi v0 = VectorXi::Zero(m);
  vector<bool> points(cord(m-1)*d,false);
  points[cord.dot(off)] = true;
  for(int j=0; j<v.size(); j++){
    VectorXi vj = v[j];
    if(vj.lpNorm<1>()>=g){  //TODO: make lower bound use the given norm
      int addr = cord.dot(vj+off);
      points[addr] = true;

      vector<VectorXi> subj = {vj};
      vector<VectorXi> valj = {v0,vj};
      futes.push_back(pool.submit(
        subsetTest, 1, subj, j+1, valj, points));
      points[addr] = false;
    }
  }

  vector<vector<VectorXi>> res = {};
  for(auto& f : futes){
    vector<vector<VectorXi>> resA = f.get();
    res.insert(res.end(), resA.begin(), resA.end());
  }
  return res;
}

int main(int argc, char** argv){
  //Take 5 command line arguments
  int a = atoi(argv[1]);    //0 -> Delta, 1 -> |A|1 norm for bounds
  int m = atoi(argv[2]);    //vector dimension
  int s = atoi(argv[3]);    //desired set size
  int g = atoi(argv[4]);    //DELTA lower bound
  int c = atoi(argv[5]);    //upper bound in norm
  vector<vector<VectorXi>> res = MaxSupp(a,m,s,g,c); //env
  std::ofstream out;
  out.open("out.txt");
  out << m << " " << s << "\n";
  for(auto S : res){
    for(auto v : S)
      for(int i=0; i<v.size(); i++)
        out << v(i) << " ";
    out << "\n";
  }
}
