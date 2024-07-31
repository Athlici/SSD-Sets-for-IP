#include <vector>
#include <fstream>
#include <iostream>
//#include <stdlib.h>
//#include "gurobi_c++.h"
#include <eigen3/Eigen/Dense>

using std::vector;
using std::to_string;
using Eigen::VectorXi;

int abs(int a){
  if(a >= 0)
    return  a;
  else
    return -a;
}

//Generates all bounded vectors of dimension m.
//For a=0 the bound c is on the largest coefficient,
//for a=1 the bound c is on the one norm of the vector.
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

int main(int argc, char** argv){
  int M = 100;
  uint8_t a = atoi(argv[1]);
  uint8_t m = atoi(argv[2]);
  uint8_t c = atoi(argv[3]);

  vector<VectorXi> v = vecs(a,m,c);
  int n = v.size()-1;

//  for(int j=0; j<m; j++){
//    for(int i=1; i<v.size(); i++)
//      std::cout << v[i](j) << " ";
//    std::cout << "\n";
//  }

  std::ofstream out;
  std::string fn = "a"+to_string(a)+"m"+to_string(m)+"c"+to_string(c);
  out.open(fn+".mat");
  out << +m << " " << n << "\n";
  for(int j=0; j<m; j++){
    for(int i=0; i<n; i++)
      out << v[i+1](j) << " ";
    out << "\n";
  }
  out.close();

  system(("graver -m "+fn).c_str());

//  std::fstream in;
//  in.open(fn+".gra");
//  int am, an;
//  in >> am;
//  in >> an;
//  VectorXi aug = VectorXi::Zero(n);
//  vector<VectorXi> augs = {};
//  bool cont = true;
//  for(int k=0; cont; k++){
//    for(int i=0; i<m; i++){
//      for(int j=0; j<n; j++)
//        cont = cont && (in >> aug(j));
//      augs.push_back(aug);
//    }
//  }

//  GRBLinExpr cl, cr;
//  double vec1[n], vecd[n], wgts[2] = {1,1};
//  char   grbI[n], grbB[n];
//  for(int i=0;i<n;i++){
//    vec1[i] = 1;
//    vecd[i] = M/10;
//    grbI[i] = GRB_INTEGER;
//    grbB[i] = GRB_BINARY;
//  }
//
//  GRBEnv *env = new GRBEnv("grb.log");
//  GRBModel *set = new GRBModel(*env);
//  set->set(GRB_IntAttr_ModelSense,GRB_MAXIMIZE);
//  GRBVar* spvars = set->addVars(NULL,vecd,vec1,grbI,NULL,n);
//  GRBVar* snvars = set->addVars(NULL,vecd,vec1,grbI,NULL,n);
//
//  for(int j=0; j<n; j++){
//    GRBVar grbv[2] = {spvars[j], snvars[j]};
//    set->addSOS(grbv,wgts,2,GRB_SOS_TYPE1);
//  }
//
//  for(int i=0; i<m; i++){
//    cl = 0; cr = 0;
//    for(int j=0; j<n; j++)
//      vecd[j] = v[j+1][i];
//    cl.addTerms(vecd,spvars,n);
//    cr.addTerms(vecd,snvars,n);
//    set->addConstr(cl, GRB_EQUAL, cr);
//  }
//
//  for(int i=0; i<am; i++){
//    GRBVar* pset = set->addVars(NULL,NULL,NULL,grbB,NULL,n);
//    GRBVar* nset = set->addVars(NULL,NULL,NULL,grbB,NULL,n);
//    for(int j=0; j<n; j++){
//      if(augs[i][j]>=0){
//        cl = spvars[j]-M*pset[j];
//        set->addConstr(cl,GRB_LESS_EQUAL,augs[i][j]-1);
//        cl = snvars[j]-M*nset[j];
//        set->addConstr(cl,GRB_LESS_EQUAL,augs[i][j]-1);
//      } else {
//        cl = snvars[j]-M*pset[j];
//        set->addConstr(cl,GRB_LESS_EQUAL,-augs[i][j]-1);
//        cl = spvars[j]-M*nset[j];
//        set->addConstr(cl,GRB_LESS_EQUAL,-augs[i][j]-1);
//      }
//    }
//    cl = 0; cr = 0;
//    cl.addTerms(vec1,pset,n);
//    cr.addTerms(vec1,nset,n);
//    set->addConstr(cl,GRB_LESS_EQUAL,n-1);
//    set->addConstr(cr,GRB_LESS_EQUAL,n-1);
//  }
//
//  set->write("test.lp");
//
//  set->optimize();
//  set->write("test.sol");
//
//
//  for(int j=0;j<n;j++){
//    int vc = std::lround(spvars[j].get(GRB_DoubleAttr_X)-snvars[j].get(GRB_DoubleAttr_X));
//    if(vc != 0){
//      std::cout << vc << "x(";
//      for(int i=0; i<m; i++)
//        std::cout << v[j+1][i] << ",";
//      std::cout << "\b) ";
//    }
//  }
//  std::cout << "\n";

}
