// rename file isr_test_io.h
#ifndef ISR_TEST_UTIL_FILE_READING_H
#define ISR_TEST_UTIL_FILE_READING_H

#include <iostream>
#include <fstream>
#include "isr_test_types.h"
#include "isr_test_normal_utils.h"
#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <boost/property_map/property_map.hpp>
//PMP
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Halfedge_index Halfedge_index;

//boost
typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;


bool is_mesh_file(const std::string input_filename)
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

bool is_point_set_file(const std::string input_filename)
{
  std::string extension = input_filename.substr(input_filename.find_last_of('.'));
  if (extension == ".xyz" || extension == ".XYZ" || 
      extension == ".pwn" || extension == ".PWN") 
  {
    return true;
  }
  return false;
}

bool read_mesh_file(const std::string &input_filename, Mesh &input_mesh) /*mieux gerer la variable success*/
{
  bool success = true; // remove
  std::ifstream stream(input_filename);
  if (!stream)
  {
    std::cerr << "Error: unable read the file " << input_filename << std::endl;
    return false;
  }

  std::string extension = input_filename.substr(input_filename.find_last_of('.'));
  if (extension == ".off" || extension == ".OFF") 
  {
    stream >> input_mesh;
    if(!input_mesh.is_valid() || input_mesh.is_empty()) // first is_empty() as it is faster, can use our is_valid()
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

bool read_point_set_file(const std::string &input_filename, PwnList &pwnl)
{
  bool success = true; // remove
  std::ifstream stream(input_filename);

  if (!stream ||
      !CGAL::read_xyz_points(
                              stream,
                              std::back_inserter(pwnl),
                              CGAL::parameters::point_map(Point_map()).
                              normal_map(Normal_map())))
    {
      std::cerr << "Error: unable to read .xyz/.pwn file" << input_filename << std::endl;
      success = false; // return false
    }

  return success; // return true
}

bool get_point_set_from_file(const std::string &in_file, PwnList &pwnl)
{
  //reads input file
  if(is_mesh_file(in_file))
  {
    Mesh input_mesh;
    if (read_mesh_file(in_file, input_mesh)) {
      Mesh::Property_map<vertex_descriptor, Vector> vnormals_pm = 
        input_mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
      compute_area_weighted_vertex_normals(input_mesh, vnormals_pm);
      BOOST_FOREACH(vertex_descriptor v, input_mesh.vertices()) {
        const Point& p = input_mesh.point(v);
        Vector n = vnormals_pm[v];
        //Vector n = PMP::compute_vertex_normal(v , input_mesh , vnormals_pm);
        pwnl.push_back(std::make_pair(p, n));
      }
    }
    else
      return false;
  }

  else if (is_point_set_file(in_file)) {
    if (!read_point_set_file(in_file, pwnl)) {
      return false;
    }
  }

  else {
    std::cerr << "Error: file not supported" << in_file << std::endl;
    return false;
  }

  //checks requirements
  std::size_t nb_points = pwnl.size();

  if (nb_points == 0) {
    std::cerr << "Error: empty point set" << std::endl;
    return false;
  }

  bool points_have_normals = (pwnl.begin()->second != CGAL::NULL_VECTOR);
  if ( ! points_have_normals ) {
    std::cerr << "Input point set not supported: this reconstruction method requires unoriented normals" << std::endl;
    return false;
  }

  return true;
}

#endif //ISR_TEST_UTIL_FILE_READING_H
