#include <vector>
#include <fstream>
#include <iostream>

using std::vector;

int supp(int a){
  return a!=0;
}

int abs(int a){
  return (a>=0) ? a : -a;
}

int (*p[2]) (int) = {supp, abs};

int main(int argc, char** argv){
  int   s  = atoi(argv[1]);
  char* fn =      argv[2] ;

  int max = 0;
  vector<vector<int>> gmax = {};

  std::fstream in;
  in.open(fn);
  int am, an;
  in >> am;
  in >> an;
  for(int i=0; i<am; i++){
    int tmp = 0;
    vector<int> g(an,0);
    for(int j=0; j<an; j++){
      in >> g[j];
      tmp += p[s](g[j]);
    }
    if(tmp > max){
      max = tmp;
      gmax.clear();
    }
    if(tmp == max)
      gmax.push_back(g);
  }

  std::cout << gmax.size() << " " << an << "\n";
  for(vector<int> & g : gmax){
    for(int j=0; j<an; j++)
      std::cout << g[j] << " ";
    std::cout << "\n";
  }
}
