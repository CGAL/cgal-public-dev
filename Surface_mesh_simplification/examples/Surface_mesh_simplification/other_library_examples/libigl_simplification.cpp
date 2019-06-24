#include <igl/readOFF.h>
#include <igl/qslim.h>
#include <igl/writeOFF.h>
#include <igl/readSTL.h>
#include <igl/remove_duplicate_vertices.h>
#include <iostream>
#include <chrono>
#include <algorithm>

int main(int argc, char* argv[]) {
  if(argc < 3) {
    std::cerr << "Usage: ./"
              << argv[0]
              << " file.off <resulting-number-of-faces>"
              << std::endl;
    return 1;
  }

  Eigen::MatrixXd temp_V, V;
  Eigen::MatrixXi F;
  if(!igl::readOFF(argv[1], temp_V, F)) {
    /*Eigen::MatrixXd N;
    Eigen::MatrixXi SVI, SVJ;
    igl::readSTL(argv[1], temp_V, F, N);
    igl::remove_duplicate_vertices(temp_V, 0, V, SVI, SVJ);
    std::for_each(F.data(),F.data()+F.size(),[&SVJ](int & f){f=SVJ(f);});
    igl::writeOFF(std::string(argv[1]) + ".off", V, F);*/
    return 0;
  };


  std::chrono::steady_clock::time_point start_time
    = std::chrono::steady_clock::now();

  Eigen::MatrixXd U;
  Eigen::MatrixXi G;
  Eigen::VectorXi J, I;

  igl::qslim(V, F, std::stoi(argv[2]), U, G, J, I);

  std::chrono::steady_clock::time_point end_time
    = std::chrono::steady_clock::now();


  std::cout << "Time elapsed: "
    << std::chrono::duration_cast<std::chrono::milliseconds>(
          end_time - start_time
        ).count() << "ms" << std::endl;


  igl::writeOFF("out.off", U, G);

  return 0;
}
