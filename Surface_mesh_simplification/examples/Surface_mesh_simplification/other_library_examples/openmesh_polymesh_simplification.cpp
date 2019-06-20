#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <iostream>
#include <chrono>

typedef OpenMesh::PolyMesh_ArrayKernelT<>              Mesh;
typedef OpenMesh::Decimater::DecimaterT<Mesh>          Decimater;
typedef OpenMesh::Decimater::ModQuadricT<Mesh>::Handle HModQuadric;

int main(int argc, char** argv) {
  if(argc < 3) {
    std::cerr
      << "Usage: ./"
      << argv[0]
      << " file.off <resulting-number-of-vertices>"
      << std::endl;
    return 1;
  }

  Mesh mesh;

  if(!OpenMesh::IO::read_mesh(mesh, argv[1])) {
    std::cerr << argv[1] << " mesh read error" << std::endl;
    return 1;
  }


  std::chrono::steady_clock::time_point start_time
    = std::chrono::steady_clock::now();

  Decimater decimater(mesh);
  HModQuadric hModQuadric;

  decimater.add(hModQuadric);
  decimater.module(hModQuadric).unset_max_err();
  decimater.initialize();
  decimater.decimate_to(std::stoi(argv[2]));

  mesh.garbage_collection();

  std::chrono::steady_clock::time_point end_time
    = std::chrono::steady_clock::now();


  std::cout << "Time elapsed: "
    << std::chrono::duration_cast<std::chrono::milliseconds>(
          end_time - start_time
        ).count() << "ms" << std::endl;

  if(!OpenMesh::IO::write_mesh(mesh, "out.off")) {
    std::cerr << "out.off mesh write error" << std::endl;
    return 1;
  }
  return 0;
}
