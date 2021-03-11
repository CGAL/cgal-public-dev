// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_CONCAVE_HULL_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_CONCAVE_HULL_2_H

#define _USE_MATH_DEFINES

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <cassert>
#include <unordered_map>
#include <cstdint>

#include <flann/flann.hpp>

#define DBL_EPSILONL 2.2204460492503131e-16

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  namespace flann_hull {

  using std::uint64_t;

  struct Point
  {
    double x = 0.0;
    double y = 0.0;
    uint64_t id = 0;

    Point() = default;

    Point(double x, double y)
      : x(x)
      , y(y)
    {}
  };

  struct PointValue
  {
    Point point;
    double distance = 0.0;
    double angle = 0.0;
  };

  static const size_t stride = 24; // size in bytes of x, y, id

  using PointVector = std::vector<Point>;
  using PointValueVector = std::vector<PointValue>;
  using LineSegment = std::pair<Point, Point>;

  // Floating point comparisons
  auto Equal(double a, double b) -> bool;
  auto Zero(double a) -> bool;
  auto LessThan(double a, double b) -> bool;
  auto LessThanOrEqual(double a, double b) -> bool;
  auto GreaterThan(double a, double b) -> bool;

  // I/O
  auto Usage() -> void;
  auto FindArgument(int argc, char **argv, const std::string &name) -> int;
  auto ParseArgument(int argc, char **argv, const std::string &name, std::string &val) -> int;
  auto ParseArgument(int argc, char **argv, const std::string &name, int &val) -> int;
  auto HasExtension(const std::string &filename, const std::string &ext) -> bool;
  auto ReadFile(const std::string &filename, int fieldX = 1, int fieldY = 2) -> PointVector;
  auto Print(const std::string &filename, const PointVector &points) -> void;
  auto Print(FILE *out, const PointVector &points, const char *format = "%.3f  %.3f\n") -> void;
  auto Split(const std::string &value, const char *delims) -> std::vector<std::string>;

  // Algorithm-specific
  auto NearestNeighboursFlann(flann::Index<flann::L2<double>> &index, const Point &p, size_t k) -> PointValueVector;
  auto ConcaveHull(PointVector &dataset, size_t k, bool iterate) -> PointVector;
  auto ConcaveHull(PointVector &dataset, size_t k, PointVector &hull) -> bool;
  auto SortByAngle(PointValueVector &values, const Point &p, double prevAngle) -> PointVector;
  auto AddPoint(PointVector &points, const Point &p) -> void;

  // General maths
  auto PointsEqual(const Point &a, const Point &b) -> bool;
  auto Angle(const Point &a, const Point &b) -> double;
  auto NormaliseAngle(double radians) -> double;
  auto PointInPolygon(const Point &p, const PointVector &list) -> bool;
  auto Intersects(const LineSegment &a, const LineSegment &b) -> bool;

  // Point list utilities
  auto FindMinYPoint(const PointVector &points) -> Point;
  auto RemoveDuplicates(PointVector &points) -> void;
  auto IdentifyPoints(PointVector &points) -> void;
  auto RemoveHull(PointVector &points, const PointVector &hull) -> PointVector::iterator;
  auto MultiplePointInPolygon(PointVector::iterator begin, PointVector::iterator end, const PointVector &hull) -> bool;

  // Unit tests
  auto TestAngle() -> void;
  auto TestIntersects() -> void;
  auto TestSplit() -> void;

  // Print program usage info.
  auto Usage() -> void
  {
    std::cout << "Usage: concave.exe filename [-out arg] [-k arg] [-field_for_x arg] [-field_for_y arg] [-no_out] [-no_iterate]\n";
    std::cout << "\n";
    std::cout << " filename      (required) : file of input coordinates, one row per point.\n";
    std::cout << " -out          (optional) : output file for the hull polygon coordinates. Default=stdout.\n";
    std::cout << " -k            (optional) : start iteration K value. Default=3.\n";
    std::cout << " -field_for_x  (optional) : 1-based column number of input for x-coordinate. Default=1.\n";
    std::cout << " -field_for_y  (optional) : 1-based column number of input for y-coordinate. Default=2.\n";
    std::cout << " -no_out       (optional) : disable output of the hull polygon coordinates.\n";
    std::cout << " -no_iterate   (optional) : stop after only one iteration of K, irrespective of result.\n";
  }

  // Get command line index of name
  auto FindArgument(int argc, char **argv, const std::string &name) -> int
  {
    for (int i = 1; i < argc; ++i)
      {
      if (std::string(argv[i]) == name)
        return i;
      }
    return -1;
  }

  // Get the command line value (string) for name
  auto ParseArgument(int argc, char **argv, const std::string &name, std::string &val) -> int
  {
    int index = FindArgument(argc, argv, name) + 1;
    if (index > 0 && index < argc)
      val = argv[index];

    return index - 1;
  }

  // Get the command line value (int) for name
  auto ParseArgument(int argc, char **argv, const std::string &name, int &val) -> int
  {
    int index = FindArgument(argc, argv, name) + 1;

    if (index > 0 && index < argc)
      val = atoi(argv[index]);

    return (index - 1);
  }

  // Check whether a string ends with a specified suffix.
  auto HasExtension(const std::string &filename, const std::string &ext) -> bool
  {
    if (filename.length() >= ext.length())
      return (0 == filename.compare(filename.length() - ext.length(), ext.length(), ext));
    return false;
  }

  // Read a file of coordinates into a vector. First two fields of comma/tab/space delimited input are used.
  auto ReadFile(const std::string &filename, int fieldX, int fieldY) -> PointVector
  {
    fieldX--; // from 1-based index to 0-based
    fieldY--;

    PointVector list;
    Point p;
    std::string line;
    std::vector<std::string> tokens;

    std::ifstream fin(filename.c_str());
    if (fin.is_open())
      {
      while (fin.good())
        {
        getline(fin, line);
        if (!line.empty())
          {
          tokens = Split(line, " ,\t");
          if (tokens.size() >= 2)
            {
            p.x = std::atof(tokens[fieldX].c_str());
            p.y = std::atof(tokens[fieldY].c_str());
            list.push_back(p);
            }
          }
        }
      }

    return list;
  }

  // Output a point list to a file
  auto Print(const std::string &filename, const PointVector &points) -> void
  {
    std::string format;
    if (HasExtension(filename, ".csv"))
      format = "%.3f,%.3f\n";
    else
      format = "%.3f  %.3f\n";

    FILE *out = fopen(filename.c_str(), "w+");
    if (out)
      {
      Print(out, points, format.c_str());

      fclose(out);
      }
  }

  // Output a point list to a stream with a given format string
  auto Print(FILE *out, const PointVector &points, const char *format) -> void
  {
    for (const auto &p : points)
      {
      fprintf(out, format, p.x, p.y);
      }
  }

  // Iteratively call the main algorithm with an increasing k until success
  auto ConcaveHull(PointVector &dataset, size_t k, bool iterate) -> PointVector
  {
    while (k < dataset.size())
      {
      PointVector hull;
      if (ConcaveHull(dataset, k, hull) || !iterate)
        {
        return hull;
        }
      k++;
      }

    return{};
  }

  // The main algorithm from the Moreira-Santos paper.
  auto ConcaveHull(PointVector &pointList, size_t k, PointVector &hull) -> bool
  {
    hull.clear();

    if (pointList.size() < 3)
    {
      return true;
    }
    if (pointList.size() == 3)
    {
      hull = pointList;
      return true;
    }

    // construct a randomized kd-tree index using 4 kd-trees
    // 2 columns, but stride = 24 bytes in width (x, y, ignoring id)
    flann::Matrix<double> matrix(&(pointList.front().x), pointList.size(), 2, stride);
    flann::Index<flann::L2<double>> flannIndex(matrix, flann::KDTreeIndexParams(4));
    flannIndex.buildIndex();

    // std::cout << "\rFinal 'k'        : " << k;

    // Initialise hull with the min-y point
    Point firstPoint = FindMinYPoint(pointList);
    AddPoint(hull, firstPoint);

    // Until the hull is of size > 3 we want to ignore the first point from nearest neighbour searches
    Point currentPoint = firstPoint;
    flannIndex.removePoint(firstPoint.id);

    double prevAngle = 0.0;
    int step = 1;

    // Iterate until we reach the start, or until there's no points left to process
    while ((!PointsEqual(currentPoint, firstPoint) || step == 1) && hull.size() != pointList.size())
      {
      if (step == 4)
        {
        // Put back the first point into the dataset and into the flann index
        firstPoint.id = pointList.size();
        flann::Matrix<double> firstPointMatrix(&firstPoint.x, 1, 2, stride);
        flannIndex.addPoints(firstPointMatrix);
        }

      PointValueVector kNearestNeighbours = NearestNeighboursFlann(flannIndex, currentPoint, k);
      PointVector cPoints = SortByAngle(kNearestNeighbours, currentPoint, prevAngle);

      bool its = true;
      size_t i = 0;

      while (its && i < cPoints.size())
        {
        size_t lastPoint = 0;
        if (PointsEqual(cPoints[i], firstPoint))
          lastPoint = 1;

        size_t j = 2;
        its = false;

        while (!its && j < hull.size() - lastPoint)
          {
          auto line1 = std::make_pair(hull[step - 1], cPoints[i]);
          auto line2 = std::make_pair(hull[step - j - 1], hull[step - j]);
          its = Intersects(line1, line2);
          j++;
          }

        if (its)
          i++;
        }

      if (its)
        return false;

      currentPoint = cPoints[i];

      AddPoint(hull, currentPoint);

      prevAngle = Angle(hull[step], hull[step - 1]);

      flannIndex.removePoint(currentPoint.id);

      step++;
      }

    // The original points less the points belonging to the hull need to be fully enclosed by the hull in order to return true.
    PointVector dataset = pointList;

    auto newEnd = RemoveHull(dataset, hull);
    bool allEnclosed = MultiplePointInPolygon(begin(dataset), newEnd, hull);

    return allEnclosed;
  }

  // Compare a and b for equality
  auto Equal(double a, double b) -> bool
  {
    return fabs(a - b) <= DBL_EPSILONL;
  }

  // Compare value to zero
  auto Zero(double a) -> bool
  {
    return fabs(a) <= DBL_EPSILONL;
  }

  // Compare for a < b
  auto LessThan(double a, double b) -> bool
  {
    return a < (b - DBL_EPSILONL);
  }

  // Compare for a <= b
  auto LessThanOrEqual(double a, double b) -> bool
  {
    return a <= (b + DBL_EPSILONL);
  }

  // Compare for a > b
  auto GreaterThan(double a, double b) -> bool
  {
    return a > (b + DBL_EPSILONL);
  }

  // Compare whether two points have the same x and y
  auto PointsEqual(const Point &a, const Point &b) -> bool
  {
    return Equal(a.x, b.x) && Equal(a.y, b.y);
  }

  // Remove duplicates in a list of points
  auto RemoveDuplicates(PointVector &points) -> void
  {
    sort(begin(points), end(points), [](const Point & a, const Point & b)
      {
      if (Equal(a.x, b.x))
        return LessThan(a.y, b.y);
      else
        return LessThan(a.x, b.x);
      });

    auto newEnd = unique(begin(points), end(points), [](const Point & a, const Point & b)
      {
      return PointsEqual(a, b);
      });

    points.erase(newEnd, end(points));
  }

  // Uniquely id the points for binary searching
  auto IdentifyPoints(PointVector &points) -> void
  {
    uint64_t id = 0;

    for (auto itr = begin(points); itr != end(points); ++itr, ++id)
      {
      itr->id = id;
      }
  }

  // Find the point having the smallest y-value
  auto FindMinYPoint(const PointVector &points) -> Point
  {
    assert(!points.empty());

    auto itr = min_element(begin(points), end(points), [](const Point & a, const Point & b)
      {
      if (Equal(a.y, b.y))
        return GreaterThan(a.x, b.x);
      else
        return LessThan(a.y, b.y);
      });

    return *itr;
  }

  // Lookup by ID and remove a point from a list of points
  auto RemovePoint(PointVector &list, const Point &p) -> void
  {
    auto itr = std::lower_bound(begin(list), end(list), p, [](const Point & a, const Point & b)
      {
      return a.id < b.id;
      });

    assert(itr != end(list) && itr->id == p.id);

    if (itr != end(list))
      list.erase(itr);
  }

  // Add a point to a list of points
  auto AddPoint(PointVector &points, const Point &p) -> void
  {
    points.push_back(p);
  }

  // Return the k-nearest points in a list of points from the given point p (uses Flann library).
  auto NearestNeighboursFlann(flann::Index<flann::L2<double>> &index, const Point &p, size_t k) -> PointValueVector
  {
    std::vector<int> vIndices(k);
    std::vector<double> vDists(k);
    double test[] = { p.x, p.y };

    flann::Matrix<double> query(test, 1, 2);
    flann::Matrix<int> mIndices(vIndices.data(), 1, static_cast<int>(vIndices.size()));
    flann::Matrix<double> mDists(vDists.data(), 1, static_cast<int>(vDists.size()));

    int count_ = index.knnSearch(query, mIndices, mDists, k, flann::SearchParams(128));
    size_t count = static_cast<size_t>(count_);

    PointValueVector result(count);

    for (size_t i = 0; i < count; ++i)
      {
      int id = vIndices[i];
      const double *point = index.getPoint(id);
      result[i].point.x = point[0];
      result[i].point.y = point[1];
      result[i].point.id = id;
      result[i].distance = vDists[i];
      }

    return result;
  }

  // Returns a list of points sorted in descending order of clockwise angle
  auto SortByAngle(PointValueVector &values, const Point &from, double prevAngle) -> PointVector
  {
    for_each(begin(values), end(values), [from, prevAngle](PointValue & to)
      {
      to.angle = NormaliseAngle(Angle(from, to.point) - prevAngle);
      });

    sort(begin(values), end(values), [](const PointValue & a, const PointValue & b)
      {
      return GreaterThan(a.angle, b.angle);
      });

    PointVector angled(values.size());

    transform(begin(values), end(values), begin(angled), [](const PointValue & pv)
      {
      return pv.point;
      });

    return angled;
  }

  // Get the angle in radians measured clockwise from +'ve x-axis
  auto Angle(const Point &a, const Point &b) -> double
  {
    double angle = -atan2(b.y - a.y, b.x - a.x);

    return NormaliseAngle(angle);
  }

  // Return angle in range: 0 <= angle < 2PI
  auto NormaliseAngle(double radians) -> double
  {
    if (radians < 0.0)
      return radians + M_PI + M_PI;
    else
      return radians;
  }

  // Return the new logical end after removing points from dataset having ids belonging to hull
  auto RemoveHull(PointVector &points, const PointVector &hull) -> PointVector::iterator
  {
    std::vector<uint64_t> ids(hull.size());

    transform(begin(hull), end(hull), begin(ids), [](const Point & p)
      {
      return p.id;
      });

    sort(begin(ids), end(ids));

    return remove_if(begin(points), end(points), [&ids](const Point & p)
      {
      return binary_search(begin(ids), end(ids), p.id);
      });
  }

  // Uses OpenMP to determine whether a condition exists in the specified range of elements. https://msdn.microsoft.com/en-us/library/ff521445.aspx
  template <class InIt, class Predicate>
  bool omp_parallel_any_of(InIt first, InIt last, const Predicate &pr)
  {
    typedef typename std::iterator_traits<InIt>::value_type item_type;

    // A flag that indicates that the condition exists.
    bool found = false;

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(last - first); ++i)
      {
      if (!found)
        {
        item_type &cur = *(first + i);

        // If the element satisfies the condition, set the flag to cancel the operation.
        if (pr(cur))
          {
          found = true;
          }
        }
      }

    return found;
  }

  // Check whether all points in a begin/end range are inside hull.
  auto MultiplePointInPolygon(PointVector::iterator begin, PointVector::iterator end, const PointVector &hull) -> bool
  {
    auto test = [&hull](const Point & p)
      {
      return !PointInPolygon(p, hull);
      };

    bool anyOutside = true;

  #if defined USE_OPENMP

    anyOutside = omp_parallel_any_of(begin, end, test); // multi-threaded

  #else

    anyOutside = std::any_of(begin, end, test); // single-threaded

  #endif

    return !anyOutside;
  }

  // Point-in-polygon test
  auto PointInPolygon(const Point &p, const PointVector &list) -> bool
  {
    if (list.size() <= 2)
      return false;

    const double &x = p.x;
    const double &y = p.y;

    int inout = 0;
    auto v0 = list.begin();
    auto v1 = v0 + 1;

    while (v1 != list.end())
      {
      if ((LessThanOrEqual(v0->y, y) && LessThan(y, v1->y)) || (LessThanOrEqual(v1->y, y) && LessThan(y, v0->y)))
        {
        if (!Zero(v1->y - v0->y))
          {
          double tdbl1 = (y - v0->y) / (v1->y - v0->y);
          double tdbl2 = v1->x - v0->x;

          if (LessThan(x, v0->x + (tdbl2 * tdbl1)))
            inout++;
          }
        }

      v0 = v1;
      v1++;
      }

    if (inout == 0)
      return false;
    else if (inout % 2 == 0)
      return false;
    else
      return true;
  }

  // Test whether two line segments intersect each other
  auto Intersects(const LineSegment &a, const LineSegment &b) -> bool
  {
    // https://www.topcoder.com/community/data-science/data-science-tutorials/geometry-concepts-line-intersection-and-its-applications/

    const double &ax1 = a.first.x;
    const double &ay1 = a.first.y;
    const double &ax2 = a.second.x;
    const double &ay2 = a.second.y;
    const double &bx1 = b.first.x;
    const double &by1 = b.first.y;
    const double &bx2 = b.second.x;
    const double &by2 = b.second.y;

    double a1 = ay2 - ay1;
    double b1 = ax1 - ax2;
    double c1 = a1 * ax1 + b1 * ay1;
    double a2 = by2 - by1;
    double b2 = bx1 - bx2;
    double c2 = a2 * bx1 + b2 * by1;
    double det = a1 * b2 - a2 * b1;

    if (Zero(det))
      {
      return false;
      }
    else
      {
      double x = (b2 * c1 - b1 * c2) / det;
      double y = (a1 * c2 - a2 * c1) / det;

      bool on_both = true;
      on_both = on_both && LessThanOrEqual(std::min(ax1, ax2), x) && LessThanOrEqual(x, std::max(ax1, ax2));
      on_both = on_both && LessThanOrEqual(std::min(ay1, ay2), y) && LessThanOrEqual(y, std::max(ay1, ay2));
      on_both = on_both && LessThanOrEqual(std::min(bx1, bx2), x) && LessThanOrEqual(x, std::max(bx1, bx2));
      on_both = on_both && LessThanOrEqual(std::min(by1, by2), y) && LessThanOrEqual(y, std::max(by1, by2));
      return on_both;
      }
  }

  // Unit test of Angle() function
  auto TestAngle() -> void
  {
    auto ToDegrees = [](double radians)
      {
      return radians * 180.0 / M_PI;
      };

    auto Test = [&](const Point &p, double expected)
      {
      double actual = ToDegrees(Angle({ 0.0, 0.0 }, p));
      assert(Equal(actual, expected));
      };

    double value = ToDegrees(atan(3.0 / 4.0));

    Test({  5.0,  0.0 }, 0.0);
    Test({  4.0,  3.0 }, 360.0 - value);
    Test({  3.0,  4.0 }, 270.0 + value);
    Test({  0.0,  5.0 }, 270.0);
    Test({ -3.0,  4.0 }, 270.0 - value);
    Test({ -4.0,  3.0 }, 180.0 + value);
    Test({ -5.0,  0.0 }, 180.0);
    Test({ -4.0, -3.0 }, 180.0 - value);
    Test({ -3.0, -4.0 }, 90.0 + value);
    Test({  0.0, -5.0 }, 90.0);
    Test({  3.0, -4.0 }, 90.0 - value);
    Test({  4.0, -3.0 }, 0.0 + value);
  }

  // Unit test the Intersects() function
  auto TestIntersects() -> void
  {
    using std::make_pair;

    std::unordered_map<char, Point> values;
    values['A'] = {  0.0,  0.0 };
    values['B'] = { -1.5,  3.0 };
    values['C'] = {  2.0,  2.0 };
    values['D'] = { -2.0,  1.0 };
    values['E'] = { -2.5,  5.0 };
    values['F'] = { -1.5,  7.0 };
    values['G'] = {  1.0,  9.0 };
    values['H'] = { -4.0,  7.0 };
    values['I'] = {  3.0, 10.0 };
    values['J'] = {  2.0, 11.0 };
    values['K'] = { -1.0, 11.0 };
    values['L'] = { -3.0, 11.0 };
    values['M'] = { -5.0,  9.5 };
    values['N'] = { -6.0,  7.5 };
    values['O'] = { -6.0,  4.0 };
    values['P'] = { -5.0,  2.0 };

    auto Test = [&values](const char a1, const char a2, const char b1, const char b2, bool expected)
      {
      assert(Intersects(make_pair(values[a1], values[a2]), make_pair(values[b1], values[b2])) == expected);
      assert(Intersects(make_pair(values[a2], values[a1]), make_pair(values[b1], values[b2])) == expected);
      assert(Intersects(make_pair(values[a1], values[a2]), make_pair(values[b2], values[b1])) == expected);
      assert(Intersects(make_pair(values[a2], values[a1]), make_pair(values[b2], values[b1])) == expected);
      };

    Test('B', 'D', 'A', 'C', false);
    Test('A', 'B', 'C', 'D', true);
    Test('L', 'K', 'H', 'F', false);
    Test('E', 'C', 'F', 'B', true);
    Test('P', 'C', 'E', 'B', false);
    Test('P', 'C', 'A', 'B', true);
    Test('O', 'E', 'C', 'F', false);
    Test('L', 'C', 'M', 'N', false);
    Test('L', 'C', 'N', 'B', false);
    Test('L', 'C', 'M', 'K', true);
    Test('L', 'C', 'G', 'I', false);
    Test('L', 'C', 'I', 'E', true);
    Test('M', 'O', 'N', 'F', true);
  }

  auto TestSplit() -> void
  {
    std::vector<double> expected = { -123.456, -987.654 };

    auto Test = [&expected](const std::string &input)
      {
      auto actual = Split(input, " ,\t");
      assert(actual.size() >= 2);
      assert(Equal(atof(actual[0].c_str()), expected[0]));
      assert(Equal(atof(actual[1].c_str()), expected[1]));
      };

    Test("-123.456 -987.654");
    Test("-123.4560 -987.6540");
    Test("-123.45600 -987.65400");
    Test("-123.456 -987.654 ");
    Test("-123.456 -987.654 100.5");
    Test("-123.456 -987.654 hello");
    Test("-123.456 -987.654 hello world");

    Test("-123.456,-987.654");
    Test("-123.4560,-987.6540");
    Test("-123.45600,-987.65400");
    Test("-123.456,-987.654,");
    Test("-123.456,-987.654,100.5");
    Test("-123.456,-987.654,hello");
    Test("-123.456,-987.654,hello,world");

    Test("-123.456\t-987.654");
    Test("-123.4560\t-987.6540");
    Test("-123.45600\t-987.65400");
    Test("-123.456\t-987.654\t");
    Test("-123.456\t-987.654\t100.5");
    Test("-123.456\t-987.654\thello");
    Test("-123.456\t-987.654\thello\tworld");

    Test(" -123.456   -987.654   ");
    Test(" -123.4560  -987.6540  ");
    Test(" -123.45600 -987.65400 ");
    Test(" -123.456   -987.654  ");
    Test(" -123.456   -987.654  100.5");
    Test(" -123.456   -987.654  hello");
    Test(" -123.456   -987.654  hello   world");
  }

  // String tokenise using any one of delimiters, adjacent spaces are treated as one
  auto Split(const std::string &value, const char *delims) -> std::vector<std::string>
  {
    std::vector<std::string> ret;

    size_t start = value.find_first_not_of(' ', 0);
    while (start != std::string::npos)
      {
      size_t pos = value.find_first_of(delims, start);
      if (pos == std::string::npos)
        {
        ret.push_back(value.substr(start));
        break;
        }
      else
        {
        ret.push_back(value.substr(start, pos - start));

        if (value[pos] == ' ' && strchr(delims, ' '))
          start = value.find_first_not_of(' ', pos);
        else
          start = pos + 1;
        }
      }

    return ret;
  }

  } // flann_hull

  template<typename GeomTraits>
  class Concave_hull_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;

    using Points_3 = std::vector<Point_3>;
    using Triangulation = internal::Triangulation<Traits>;
    using Partition_traits = CGAL::Partition_traits_2<Traits>;
    using Polygon_2 = typename Partition_traits::Polygon_2;
    using Location_type = typename Triangulation::Delaunay::Locate_type;

    Concave_hull_2(const Points_3& points) :
    m_points(points)
    { }

    void create(
      Points_3& output) {

      flann_hull::PointVector input;
      input.reserve(m_points.size());
      for (const auto& p : m_points)
        input.push_back(flann_hull::Point(p.x(), p.y()));

      bool iterate = true; int k = 3;
      flann_hull::RemoveDuplicates(input);
      flann_hull::IdentifyPoints(input);

      std::cout << "Hull input: " << input.size() << std::endl;

      flann_hull::PointVector hull =
        flann_hull::ConcaveHull(input, (size_t) k, iterate);

      std::cout << "Hull output: " << hull.size() << std::endl;

      output.clear();
      output.reserve(hull.size());
      for (const auto& p : hull)
        output.push_back(Point_3(p.x, p.y, FT(0)));
    }

    template<typename Pixel>
    void set_interior_labels(
      const Points_3& hull,
      std::vector<Pixel>& point_cloud) {

      Polygon_2 polygon;
      for (std::size_t i = 0; i < hull.size() - 1; ++i) {
        const auto& p = hull[i];
        polygon.push_back(Point_2(p.x(), p.y()));
      }

      std::vector<Polygon_2> convex;
      CGAL::approx_convex_partition_2(
        polygon.vertices_begin(), polygon.vertices_end(),
        std::back_inserter(convex));

      std::vector<Triangulation> tris;
      tris.resize(convex.size());

      for (std::size_t i = 0; i < convex.size(); ++i) {
        const auto& poly = convex[i];

        tris[i].delaunay.clear();
        for (std::size_t j = 0; j < poly.size(); ++j)
          tris[i].delaunay.insert(poly[j]);
      }

      Location_type type; int stub;
      for (auto& pixel : point_cloud) {
        const Point_2 p = Point_2(pixel.point.x(), pixel.point.y());

        bool found = false;
        for (const auto& tri : tris) {
          const auto fh = tri.delaunay.locate(p, type, stub);

          if (!tri.delaunay.is_infinite(fh)) {
            pixel.is_interior = true;
            found = true; break;
          }
        }

        if (!found)
          pixel.is_interior = false;
      }
    }

  private:
    const Points_3& m_points;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_CONCAVE_HULL_2_H
