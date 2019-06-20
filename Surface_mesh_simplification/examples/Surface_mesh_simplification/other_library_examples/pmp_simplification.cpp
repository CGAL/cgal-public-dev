#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/SurfaceSimplification.h>
#include <iostream>
#include <chrono>

int main(int argc, char *argv[]) {
  if(argc < 3) {
    std::cerr
      << "Usage: ./"
      << argv[0]
      << " file.off <resulting-number-of-vertices>"
      << std::endl;
    return 1;
  }

  pmp::SurfaceMesh mesh;
  mesh.read(argv[1]);


  std::chrono::steady_clock::time_point start_time
    = std::chrono::steady_clock::now();

  pmp::SurfaceSimplification simplification(mesh);
  simplification.simplify(std::stoi(argv[2]));

  std::chrono::steady_clock::time_point end_time
    = std::chrono::steady_clock::now();


  std::cout << "Time elapsed: "
    << std::chrono::duration_cast<std::chrono::milliseconds>(
          end_time - start_time
        ).count() << "ms" << std::endl;

  mesh.write("out.off");
  return 0;
}
