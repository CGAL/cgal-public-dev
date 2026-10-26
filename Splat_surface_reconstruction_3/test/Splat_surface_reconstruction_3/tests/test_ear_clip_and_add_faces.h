#pragma once
#include "common.h"

inline bool test_ear_clip_and_add_faces_single_tri() {
  std::vector<Point> points = {
    Point(1.0, 1.0, 0.0),
    Point(-1.0, -1.0, 0.0),
    Point(1.0, -1.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh = ear_clip_and_add_faces_helper(points, normals);

  assert(num_vertices(mesh) == 3);
  assert(num_edges(mesh) == 3);
  assert(mesh.number_of_faces() == 1);

  for (auto f : faces(mesh)) {
    std::size_t face_degree = 0;
    for (auto h : halfedges_around_face(halfedge(f, mesh), mesh))
      ++face_degree;

    assert(face_degree == 3);
  }

  assert(CGAL::is_valid_polygon_mesh(mesh));

  return true;
}

inline bool test_ear_clip_and_add_faces_square() {
  std::vector<Point> points = {
    Point(-1.0, 1.0, 0.0),
    Point(-1.0, -1.0, 0.0),
    Point(1.0, -1.0, 0.0),
    Point(1.0, 1.0, 0.0)
  };
  std::vector<Vector_3> normals = {
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0),
    Vector_3(0.0, 0.0, 1.0)
  };

  Polyhedron mesh = ear_clip_and_add_faces_helper(points, normals);

  std::cout << mesh << std::endl;

  assert(num_vertices(mesh) == 4);
  assert(num_edges(mesh) == 5);
  assert(mesh.number_of_faces() == 2);

  for (auto f : faces(mesh)) {
    std::size_t face_degree = 0;
    for (auto h : halfedges_around_face(halfedge(f, mesh), mesh))
      ++face_degree;

    assert(face_degree == 3);
  }

  assert(CGAL::is_valid_polygon_mesh(mesh));

  return true;
}


inline bool test_ear_clip_and_add_faces_wrapping_shape() {

  // This test checks the ear clipping and face addition for a shape that wraps around, 
  // and thus has opposing normals
  double n = -1.0;
  double p = 1.0;
  std::vector<Point> points = {
    Point(p,p,n),
    Point(p,p,p),
    Point(n,p,p),
    Point(n,n,p),
    Point(p,n,p),
    Point(p,n,n),
    Point(n,n,n),
    Point(n,p,n),
  };
  std::vector<Vector_3> normals = {
    normalize(Vector_3(1.0, 0.0, -1.0)),
    normalize(Vector_3(1.0, 0.0, 1.0)),
    normalize(Vector_3(0.0, 0.0, 1.0)),
    normalize(Vector_3(0.0, 0.0, 1.0)),
    normalize(Vector_3(1.0, 0.0, 1.0)),
    normalize(Vector_3(1.0, 0.0, -1.0)),
    normalize(Vector_3(0.0, 0.0, -1.0)),
    normalize(Vector_3(0.0, 0.0, -1.0))
  };

  Polyhedron mesh = ear_clip_and_add_faces_helper(points, normals);

  std::cout << mesh << std::endl;

  assert(num_vertices(mesh) == 4);
  assert(num_edges(mesh) == 5);
  assert(mesh.number_of_faces() == 2);

  for (auto f : faces(mesh)) {
    std::size_t face_degree = 0;
    for (auto h : halfedges_around_face(halfedge(f, mesh), mesh))
      ++face_degree;

    assert(face_degree == 3);
  }

  assert(CGAL::is_valid_polygon_mesh(mesh));

  return true;
}


inline void run_ear_clip_and_add_faces_tests(TestData& test_data) {
  test(test_ear_clip_and_add_faces_single_tri, "test_ear_clip_and_add_faces_single_tri", test_data);
  test(test_ear_clip_and_add_faces_square, "test_ear_clip_and_add_faces_square", test_data);
  // test_ear_clip_and_add_faces_wrapping_shape currently fails (was not run before the split)
}
