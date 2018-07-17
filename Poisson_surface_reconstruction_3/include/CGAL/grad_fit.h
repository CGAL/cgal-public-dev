#ifndef GRAD_FIT_H
#define GRAD_FIT_H
#include <iostream>
#include <set>
#include <cmath>
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Vector_3.h>
#include<CGAL/Point_3.h>
#include "random.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_3<K> Point;
typedef CGAL::Vector_3<K> Vector;

template <typename Vertex_handle, typename Tr>
Vector grad_fit(Tr* tr, std::set<Vertex_handle> vertices, Vertex_handle v)
{
  int m = vertices.size() + 1;

  //TODO: if (m < 10)
  std::vector<Point> points;//(m);
  std::vector<double> function_values;//(m);

  Point query = v->point();
  points.push_back(v->point());
  function_values.push_back(v->f());

  int i = 1;
  for(typename std::set<Vertex_handle>::iterator it = vertices.begin();
    it != vertices.end(); it++, i++)
  {
    if(tr->is_infinite(*it))
      continue;
    Vertex_handle v = *it;
    points.push_back(v->point());
    function_values.push_back(v->f());
  }

  m = points.size();
  double x, y, z, f;

  while(m < 10){
    Vector v = random_barycentric_coords<Vector>();
    int i = std::floor(random_double(0, m - 1));
    x = v[0] * points[i][0] + v[1] * points[(i + 1) % m][0]  + v[2] * points[(i + 2) % m][0];
    y = v[0] * points[i][1] + v[1] * points[(i + 1) % m][1]  + v[2] * points[(i + 2) % m][1];
    z = v[0] * points[i][2] + v[1] * points[(i + 1) % m][2]  + v[2] * points[(i + 2) % m][2];
    Point p(x, y, z);
    f = v[0] * function_values[i] + v[1] * function_values[(i + 1) % m] +
    v[2] * function_values[(i + 2) % m];
    points.push_back(p);
    function_values.push_back(f);

    m++;
  }

  double scaling_factor = -100;

//translating such that the first point is the origin
  for(std::vector<Point>::iterator it = points.begin(); it != points.end();
    it++)
  {

    x = (*it)[0] - query[0];
    y = (*it)[1] - query[1];
    z = (*it)[2] - query[2];
    *it = Point(x,y,z);

    scaling_factor = std::max(scaling_factor, std::max(std::abs(x),
    std::max(std::abs(y), std::abs(z))));
  }

  //scaling all the coordinates:
  for(std::vector<Point>::iterator it = points.begin(); it != points.end();
  it++)
  {
    x = (*it)[0] / scaling_factor;
    y = (*it)[1] / scaling_factor;
    z = (*it)[2] / scaling_factor;
    *it = Point(x,y,z);
  }
  double alpha = 15;

//Weight matrix
  Eigen::MatrixXd W(m, m);
  W.setZero();

  for(int i = 0; i < points.size(); i++)
  {
    Vector v(points[i], query);
    double dist = std::sqrt(v * v);
    W(i, i) = std::exp(-alpha * dist);
  }

  Eigen::MatrixXd A(m, 10);
  Eigen::VectorXd b(m); //function values (RHS of the equation)
  A(0, 0) = 0.0;
  A(0, 1) = 0.0;
  A(0, 2) = 0.0;
  A(0, 3) = 0.0;
  A(0, 4) = 0.0;
  A(0, 5) = 0.0;
  A(0, 6) = 0.0;
  A(0, 7) = 0.0;
  A(0, 8) = 0.0;
  A(0, 9) = 1.0;
  b(0) = function_values[0];
  for(int i = 1; i < m; i++)
  {
    A(i, 0) = points[i][0] * points[i][0];
    A(i, 1) = points[i][1] * points[i][1];
    A(i, 2) = points[i][2] * points[i][2];
    A(i, 3) = points[i][0] * points[i][1];
    A(i, 4) = points[i][1] * points[i][2];
    A(i, 5) = points[i][0] * points[i][2];
    A(i, 6) = points[i][0];
    A(i, 7) = points[i][1];
    A(i, 8) = points[i][2];
    A(i, 9) = 1.0;
    b(i) = function_values[i];
  }


//adding the weighting
  A = W * A;
  b = W * b;

  Eigen::VectorXd M = A.colPivHouseholderQr().solve(b);


  x = M(6)/scaling_factor;
  y = M(7)/scaling_factor;
  z = M(8)/scaling_factor;
  return Vector(x, y, z);


}

#endif