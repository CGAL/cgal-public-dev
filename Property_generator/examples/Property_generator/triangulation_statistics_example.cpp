#include <vector>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <CGAL/iterator.h>
#include <CGAL/algorithm.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/properties/triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_2                                      Point;
typedef CGAL::Delaunay_triangulation_2<Kernel>               Delaunay;
typedef CGAL::Creator_uniform_2<double, Point>               Creator;


template <typename Iterator>
class No_deref_iterator
#ifndef DOXYGEN_RUNNING
  : public boost::iterator_adaptor<
        No_deref_iterator<Iterator>    // Derived
      , Iterator                       // Base
      , Iterator                       // Value
      , boost::forward_traversal_tag   // Traversal type
      , Iterator                       // Reference
    >
#endif
{
 private:
    struct enabler {};

 public:
    No_deref_iterator()
        : No_deref_iterator::iterator_adaptor_() {}

    explicit No_deref_iterator(Iterator it)
        : No_deref_iterator::iterator_adaptor_(it) {}

 private:
    friend class boost::iterator_core_access;
 
    typename No_deref_iterator::reference
    dereference()
    const
    {
        return this->base();
    }
};

template<typename T>
inline No_deref_iterator<T> make_no_deref_iterator(T iterator){
    return No_deref_iterator<T>(iterator);
}

template <typename Iterator, typename Function>
double mean_result(Iterator begin, Iterator end, Function f)
{
  typename std::iterator_traits<Iterator>::difference_type count = 0;
  double sum = 0.;

  for (; begin != end; ++begin, ++count)
    sum += CGAL::to_double(f( *begin ));

  return sum / double(count);
}

template <typename Iterator, typename Function>
typename CGAL::cpp11::result_of<
    Function(typename std::iterator_traits<Iterator>::value_type)
>::type
max_result(Iterator begin, Iterator end, Function f)
{
  typedef typename CGAL::cpp11::result_of<Function(
      typename std::iterator_traits<Iterator>::value_type)>::type result_type;

  // Set max to value of first element.
  result_type max = f(*begin);

  for (++begin; begin != end; ++begin)
  {
    result_type value = f(*begin);
    if (value > max)
      max = value;
  }

  return max;
}




// Alias the namespaces to save typing.
namespace Properties = CGAL::Properties::Triangulation_2;

int main()
{
    // Number of points to generate.
    unsigned           n = 1<<10;
    Delaunay           dt;
    std::vector<Point> points;

    // Random triangulation.
    CGAL::Random_points_in_square_2<Point,Creator> g(0.5);
    CGAL::cpp11::copy_n(g, n, std::back_inserter(points));
    dt.insert(points.begin(), points.end());

    // The property functions take handle types, so we wrap the iterators.
    No_deref_iterator<Delaunay::Finite_faces_iterator> begin, end;

    begin = make_no_deref_iterator(dt.finite_faces_begin());
    end   = make_no_deref_iterator(dt.finite_faces_end());

    // Initialise the functors that we shall use.
    Properties::Area<Delaunay>         area(dt);
    Properties::Circumradius<Delaunay> circumradius(dt);
    Properties::Aspect_ratio<Delaunay> aspect_ratio(dt);
    Properties::Max_angle<Delaunay>    max_angle(dt);

    // Display some statistics about the triangulation.
    std::cout
        << "-- Information about the triangulation --" 
        << std::endl << std::left << std::setw(50)
        
        << "Mean face area"
        << mean_result(begin, end, area)
        << std::endl << std::left << std::setw(50)

        << "Largest face area"
        << max_result(begin, end, area)
        << std::endl << std::left << std::setw(50)

        << "Mean circumradius"
        << mean_result(begin, end, circumradius)
        << std::endl << std::left << std::setw(50)

        << "Maximum aspect ratio"
        << max_result(begin, end, aspect_ratio)
        << std::endl << std::left << std::setw(50)

      // << "Correlation between area and aspect ratio"
      //  << CGAL::pearson(begin, end, area, aspect_ratio)
      //  << std::endl << std::left << std::setw(50)

      //  << "Number of angles larger than three"
      //  <<  CGAL::count_result_in_interval(begin, end, max_angle, 3,4)
        << std::endl;
}


