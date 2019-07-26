#ifndef ISR_TEST_UTIL_FILE_READING_H
#define ISR_TEST_UTIL_FILE_READING_H

#include <iostream>
#include <fstream>
#include "isr_test_types.h"
#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/OBJ_reader.h>

#include <CGAL/IO/read_xyz_points.h>
#include <boost/property_map/property_map.hpp>


typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Halfedge_index Halfedge_index;



bool is_mesh(const std::string input_filename)
{
  std::string extension = input_filename.substr(input_filename.find_last_of('.'));
  if (extension == ".off" || extension == ".OFF" || 
      extension == ".stl" || extension == ".STL" || 
      extension == ".obj" || extension == ".OBJ") 
  {
    return true;
  }
  return false;
}

bool is_point_set(const std::string input_filename)
{
  std::string extension = input_filename.substr(input_filename.find_last_of('.'));
  if (extension == ".xyz" || extension == ".XYZ" || 
      extension == ".pwn" || extension == ".PWN") 
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
      return success;
    }
  }

  else
  {
    std::vector<Point> points;
    std::vector<std::vector<size_t>> faces;

    if((extension == ".stl" || extension == ".STL")  && !CGAL::read_STL(stream, points, faces, false))
    {
        std::cerr << "Error : unable to read the stl file." << std::endl;   
        success = false;
        return success;
    }
    else if((extension == ".obj" || extension == ".OBJ")  && !CGAL::read_OBJ(stream, points, faces))
    {
       std::cerr << "Error : unable to read the obj file." << input_filename << std::endl;   
       success = false;
       return success;
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
 return success; 
}

bool read_input_point_set_file(const std::string &input_filename, PwnList &pwnl)
{
  bool success = true;
  std::ifstream stream(input_filename);

  if (!stream ||
      !CGAL::read_xyz_points(
                              stream,
                              std::back_inserter(pwnl),
                              CGAL::parameters::point_map(Point_map()).
                              normal_map(Normal_map())))
    {
      std::cerr << "Error: unable to read .xyz/.pwn file" << input_filename << std::endl;
      success = false;
    }

  return success;
}
#endif //ISR_TEST_UTIL_FILE_READING_H