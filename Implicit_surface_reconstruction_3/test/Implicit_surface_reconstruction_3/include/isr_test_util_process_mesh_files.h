#ifndef ISR_TEST_UTIL_PROCESS_MF_H
#define ISR_TEST_UTIL_PROCESS_MF_H

#include <iostream>
#include <fstream>
#include "isr_test_types.h"
#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/OBJ_reader.h> 

typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Halfedge_index Halfedge_index;

bool is_mesh(std::string input_filename)
{
  std::string extension = input_filename.substr(input_filename.find_last_of('.'));
  if (extension == ".off" || extension == ".OFF" || extension == ".stl" || extension == ".STL" || extension == ".obj" || extension == ".OBJ") 
  {
    return true;
  }
  return false;
}

bool read_input_mesh_file(const std::string &input_filename, Mesh &input_mesh) /*mieux gerer la variable success*/
{
  bool success = true;
  std::ifstream stream(input_filename);
  std::string extension = input_filename.substr(input_filename.find_last_of('.'));

  if (extension == ".off" || extension == ".OFF") 
  {
    stream >> input_mesh;
    if(!stream || !input_mesh.is_valid() || input_mesh.is_empty())
    {
      std::cerr << "Error: unable read the off file " << input_filename << std::endl;
      success = false;
      return (success);
    }
  }

  else
  {
    std::vector<Point> points;
    std::vector<std::vector<size_t>> faces;

    if((extension == ".stl" || extension == ".STL")  && !CGAL::read_STL(stream, points, faces, false))
    {
        std::cerr << "Error : unable to read the stl file." << std::endl;   
       return false;
    }
    else if((extension == ".obj" || extension == ".OBJ")  && !CGAL::read_OBJ(stream, points, faces))
    {
       std::cerr << "Error : unable to read the obj file." << std::endl;   
       return false;   
    }

    for(const Point &pt : points)
    {
        input_mesh.add_vertex(pt);
    }
    for(const std::vector<size_t> &face : faces)
    {
       std::list<Vertex_index> vertices;
       for(size_t vi : face)
       {
         vertices.push_back(static_cast<Vertex_index>(vi));
       }
       input_mesh.add_face(vertices);
     }
  } 
 return true; 
}

#endif //ISR_TEST_UTIL_PROCESS_MF_H