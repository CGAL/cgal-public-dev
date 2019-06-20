#include <igl/readOFF.h>
#include <igl/qslim.h>
#include <igl/writeOFF.h>
#include <iostream>
#include <chrono>

int main(int argc, char* argv[]) {
  if(argc < 3) {
    std::cerr << "Usage: ./"
              << argv[0]
              << " <file.off> resulting-number-of-faces"
              << std::endl;
    return 1;
  }

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOFF(argv[1], V, F);

  //igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F);
  //viewer.launch();
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
