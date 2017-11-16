#include <iostream>

#include "Test_algo1test.cpp"
using namespace std;

int main(int argc, char** argv)
{
  //cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

  int geometric_weights=1;
  int boundary=0;

  if (argv[2])
    geometric_weights=atoi(argv[2]);
  if (argv[3])
    boundary=atoi(argv[3]);
  Test_algo1test bt;
  if (argv[1])
    bt.test(argv[1],geometric_weights,boundary);
  else
    bt.test("sphere2L.off",geometric_weights,boundary);
  return 0;
}
