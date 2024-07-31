#include <vector>
#include <future>
#include <fstream>
#include <iostream>
#include "gurobi_c++.h"
#include <eigen3/Eigen/Dense>
#include "BS_thread_pool.hpp"

using std::vector;
using std::future;
using std::string;
using std::to_string;
using Eigen::MatrixXi;
using Eigen::VectorXi;

//Try to find an objective such that 1 the k-th SSD set A forms a vertex.
VectorXi findObj(GRBEnv *env, int k, MatrixXi A){
  int n = A.rows(), m = A.cols();
  VectorXi sgn = VectorXi::Ones(n);
  vector<VectorXi> augs = {};       //store known augmentation vectors
  GRBVar grbv[2];                   //Helpers for initializing GuRoBi
  double vec1[n], wgts[2] = {1,1};
  char   grbI[n], grbB[n];
  for(int i=0;i<n;i++){
    vec1[i] = 1;
    grbI[i] = GRB_INTEGER;
    grbB[i] = GRB_BINARY;
  }

  //Increments signs until most significant bit is reached.
  //Can end then because of point symmetry.
  while(sgn(0)==1){
//std::cout << "Testing signs: ";
//for(int i=0;i<n;i++)
//  std::cout << sgn[i] << " ";
//std::cout << "\n ";
    GRBModel *obj = new GRBModel(*env), *aug = new GRBModel(*env);
    aug->set(GRB_IntParam_Threads, 1);      //More threads reduced overall efficiency
    aug->set(GRB_DoubleParam_MIPGap, 0.5);  //Approximate optimality is good enough
    obj->set(GRB_IntParam_Threads, 1);
    GRBVar* opvars = obj->addVars(NULL,NULL,vec1,NULL,NULL,n); //Objective positive part
    GRBVar* onvars = obj->addVars(NULL,NULL,vec1,NULL,NULL,n); //Objective negative part
    GRBVar* apvars = aug->addVars(NULL,NULL,vec1,grbI,NULL,n); //Augmentation positive part
    GRBVar* anvars = aug->addVars(NULL,vec1,vec1,grbB,NULL,n); //Augmentation negative part

    for(int i=0; i<n; i++){ //Prevent positive and negative part in augmentation variable
      grbv[0] = apvars[i];
      grbv[1] = anvars[i];
      aug->addSOS(grbv,wgts,2,GRB_SOS_TYPE1);
    }
    GRBLinExpr cl = 0, cr;
    cl.addTerms(vec1,apvars,n);
    cl.addTerms(vec1,anvars,n);
    aug->addConstr(cl, GRB_GREATER_EQUAL, 1); //At least one non-zero entry in augmentation
    double vecd[n];
    for(int j=0; j<m; j++){
      VectorXi con = sgn.cwiseProduct(A.col(j));
      for(int i=0; i<n; i++)
        vecd[i] = con[i];
      cl = 0; cr = 0;
      cl.addTerms(vecd,apvars,n);
      cr.addTerms(vecd,anvars,n);
      aug->addConstr(cl, GRB_EQUAL, cr);      //A.aug = 0
    }

    for(VectorXi imp : augs){
      VectorXi simp = sgn.cwiseProduct(imp);
      if(simp.minCoeff()>=-1){
        cl = 0; cr = -1;
        for(int i=0; i<n; i++)
          vecd[i] = simp[i];
        cl.addTerms(vecd,opvars,n);
        cr.addTerms(vecd,onvars,n);
        obj->addConstr(cl, GRB_LESS_EQUAL, cr); //aug.c<=-1 for all known augmentations
      }
    }

    //Repeat until one of the models becomes infeasible
    while(true){
      obj->optimize();  //Try to find an objective function
      int status = obj->get(GRB_IntAttr_Status);
      if(status!=GRB_OPTIMAL)
        break;

      cl = 0;
      //VectorXi dir = VectorXi::Zero(n);
      Eigen::VectorXd dir = Eigen::VectorXd::Zero(n);
      for(int i=0; i<n; i++){   //Extract found objective function
        dir(i) = opvars[i].get(GRB_DoubleAttr_X)-onvars[i].get(GRB_DoubleAttr_X);
        cl += dir(i)*apvars[i]-dir(i)*anvars[i];
      }
//std::cout << "Setting objective: ";
//for(int i=0;i<n;i++)
//  std::cout << dir[i] << " ";
//std::cout << "\n ";

      //Add and later delete the objective to the augmentation model
      auto dircon = aug->addConstr(cl, GRB_GREATER_EQUAL, 0);

      aug->optimize();  //Find an improving augmentation
      status = aug->get(GRB_IntAttr_Status);
      if(status!=GRB_OPTIMAL){  //If the is no augmentation output the sign vector
        std::cout << k << ": {";
        for(int i=0; i<n; i++)
          std::cout << sgn(i) << ",";
        std::cout << "\b} -> {";
        for(int i=0; i<n; i++){ //And find a short integer objective function
          opvars[i].set(GRB_CharAttr_VType, GRB_INTEGER); 
          onvars[i].set(GRB_CharAttr_VType, GRB_INTEGER); 
        }
        obj->optimize();    //TODO: Might want to demand optimality here?
        for(int i=0; i<n; i++)  //Output the objective function as well
          std::cout << opvars[i].get(GRB_DoubleAttr_X)-onvars[i].get(GRB_DoubleAttr_X) << ",";
        std::cout << "\b}\n";
        return sgn;
      }

      cl = 0;
      VectorXi imp = VectorXi::Zero(n);
      for(int i=0; i<n; i++){
        imp(i) = std::lround(apvars[i].get(GRB_DoubleAttr_X)-anvars[i].get(GRB_DoubleAttr_X));
        cl += imp(i)*opvars[i]-imp(i)*onvars[i];
      }
//std::cout << "Adding augmentation: ";
//for(int i=0;i<n;i++)
//  std::cout << imp[i] << " ";
//std::cout << "\n ";
      augs.push_back(sgn.cwiseProduct(imp));  //Otherwise add augmentation to known ones

      obj->addConstr(cl, GRB_LESS_EQUAL, -1);
      aug->remove(dircon);
    }

    for(int i=n-1; i>=0; i--){  //Binary style incrementing of sign vector
      sgn(i) *= -1;
      if(sgn(i) == -1)
        break;
    }

    delete [] opvars; delete [] apvars; //Plug memory holes
    delete [] onvars; delete [] anvars; //TODO: What happens on early exit?
    delete obj; delete aug;
  }

  return VectorXi::Zero(n); //Failure case
}

int main(int argc, char** argv){
  GRBEnv *env = new GRBEnv("grb.log");
  env->set(GRB_IntParam_OutputFlag, 0);
//  env->set(GRB_IntParam_LogToConsole, 0);

  BS::thread_pool pool(6);
  vector<future<VectorXi>> futes = {};

  std::fstream insets;
  insets.open(argv[1]); //Input
  int n, m, a;
  insets >> m;
  insets >> n;
  MatrixXi A = MatrixXi::Zero(n,m);
  vector<MatrixXi> matr = {};
  bool cont = true;
  for(int k=0; cont; k++){  //Read matrices from file
    for(int i=0; i<n; i++)
      for(int j=0; j<m; j++){
        cont = cont && (insets >> a);
        A(i,j) = a;
      }
    if(cont){
      matr.push_back(A);
      futes.push_back(pool.submit(findObj, env, k, A));
    }
  }
  for(int i=0; i<futes.size(); i++){  //Primitive progress tracker
    futes[i].get();
    if(i%100==0)
      std::cout << "Finished testing till: " << i << "\n";
  }
}
