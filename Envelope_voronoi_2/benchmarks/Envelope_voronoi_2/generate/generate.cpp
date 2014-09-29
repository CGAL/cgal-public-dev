/*
 * This program is used to generate different points sets.
 */
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/random_selection.h>
#include <CGAL/rational_rotation.h>
#include <CGAL/Gmpq.h>

#include <boost/program_options.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <vector>
#include <algorithm>

using namespace CGAL;
namespace po = boost::program_options;

typedef Simple_cartesian<int>            R;
typedef R::Point_2                       Point;
typedef Creator_uniform_2<double,Point>  Creator;
typedef std::vector<Point>               Points_vector;

typedef Simple_cartesian<CGAL::Gmpq>     Exact_k;
typedef Exact_k::FT FT;
typedef Exact_k::Point_2                 Exact_point_2;
typedef std::vector<Exact_point_2>       Exact_point_2_vector;

typedef Exact_k::Point_3                 Exact_point_3;
typedef std::vector<Exact_point_3>       Exact_point_3_vector;

typedef Exact_k::Plane_3                 Exact_plane_3;
typedef std::vector<Exact_plane_3>       Exact_plane_3_vector;

struct Less_then_planes
{
  bool operator () (const Exact_plane_3& a, const Exact_plane_3& b)
    {
      typedef Exact_k::Vector_3                   Vector_3;
      typedef Exact_k::Point_3                    Point_3;
      Vector_3 v1 = CGAL::orthogonal_vector(a);
      Vector_3 v2 = CGAL::orthogonal_vector(b);

      Point_3 p1(v1.x(), v1.y(), v1.z());
      Point_3 p2(v2.x(), v2.y(), v2.z());

      return CGAL::lexicographically_xyz_smaller(p1, p2);
    }
};

// The function receive an approximate angle (tau) on a circle, and
// returns rational cos and sin such that they are exact cos and sin
// of an angle close to tau.
// THIS FUNCTION SHOULD BE REPLACED BY CGAL::rational_rotation_approximation
void get_exact_angles_on_circle(double tau, FT & cos, FT & sin,
                                long round = 100)
{
  // 3) Compute a rational approximation tau' for tau.
  // 4) Compute cos_tau' = (1 - (tau')^2) / (1 + (tau')^2)
  // 5) Compute sin_tau' = 2*tau' / (1 + (tau')^2)
  // 6) Now x' = r * cos_tau', y' = r * sin_tau'
  FT taut(long(tan(tau/2) * round), round);
  FT taut_square = CGAL::square(taut);
  cos = (1 - taut_square) / (1 + taut_square);
  sin = 2*taut / (1 + taut_square);
}

void copy_to_exact(Points_vector &from, Exact_point_2_vector &to)
{
  Points_vector::iterator it = from.begin();
  for (; it != from.end(); ++it)
  {
    to.push_back(Exact_point_2(it->x(), it->y()));
  }
}

int main(int argv, char **argc)
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("points,p", po::value<int>()->default_value(100),
     "specify the number of points to generate")
    ("type,t", po::value<char>()->default_value('r'),
     "specify the type of set to produce:\n"
     "Point-set:\n"
     "  [r]andom: \tgenerate random point set in rectangle\n"
     "  [d]anis example: \tgenerate random point set on y = |x|\n"
     "  [l]ine: \tgenerate pointa on the x-axis\n"
     "  [c]ircle: \tgenerate random points on circle with 10 radius\n"
     "  [o]n circle: \tgenerate points on circle with 10 radius\n\n"

     "Circle set:\n"
     "  [R]andom: \tgenerate random circle set in rectangle\n\n"
     "  [A]nnulus: \tgenerate circles for annulus computation."
     " The centers of the circles are pertubed points on a circle.\n"
     "  [D]egenerate annulus: \tgenerate circles for degenerate annulus "
        "computation.\n\n"

     "On sphere:\n"
     "  [S]phere random: \tgenerate random points on the unit sphere\n"
     "  [P]lanes random: \tgenerate planes to intersect the unit sphere "
     "(for power diagram on the sphere, and not uniform.)\n"
     "  d[E]generate planes: \tgenerate planes to intersect the unit sphere "
     " so that all planes intersect in one point (for power diagram on the "
     "sphere, and not uniform.)\n\n"
      )
    ("sort,s", "sort the points lexicographically")
    ("shuffle", "random shuffle the points lexicographically")
    ;

  po::positional_options_description p;
  p.add("points", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argv, argc).options(desc).
            positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << "Generate a set of points in a desired configuration" <<
      std::endl;
    std::cout << desc << "\n";
    return 1;
  }

  char type;
  if (vm.count("type"))
  {
    type = vm["type"].as<char>();
  }
  else
  {
    std::cerr << "Error! No type." << std::endl;
    return -1;
  }

  std::size_t n_points;
  if (vm.count("points"))
  {
    n_points = vm["points"].as<int>();
  }
  else
  {
    std::cerr << "Error! No type." << std::endl;
    return -1;
  }

  // Create test point set. Prepare a vector for n points.
  Points_vector points;
  std::vector<int> weights;
  points.reserve(n_points);

  Exact_point_2_vector out_points;
  Exact_point_3_vector out_points_on_sphere;
  Exact_plane_3_vector out_planes_3;

  Random_points_in_square_2<Point, Creator> rnd_in_sq(n_points);
  Random_points_on_segment_2<Point, Creator> rnd_on_seg
    (Point(-n_points * 5, 0), Point(n_points * 5, 0));

  // for circle points
  Points_vector::iterator it;

  switch (type)
  {
  case 'r': // random in square
    CGAL::copy_n(rnd_in_sq, n_points, std::back_inserter(points));
    copy_to_exact(points, out_points);
    break;

  case 'd': // dani's example
  {
    std::size_t i = 0;

    if (n_points % 2 == 1)
    {
      points.push_back(Point(0, 0));
      i = 1;
    }
    for (; i*2 < n_points; ++i)
    {
      points.push_back(Point(i+1, i+1));
      points.push_back(Point(-i-1, i+1));
    }

    copy_to_exact(points, out_points);
    break;
  }
  case 'l': // creates points on the x-axis
  {
    for (std::size_t i = 0; i < n_points; ++i)
    {
      points.push_back(Point(i, 0));
    }
    copy_to_exact(points, out_points);
    break;
  }
  case 'c': // random on circle with radius 10
  {
    Random_points_on_circle_2<Point, Creator> rnd_on_circ(1000);
    CGAL::copy_n(rnd_on_circ, n_points, std::back_inserter(points));

    // now we have points who are not exactly on the sphere. We'll use
    // them as directions, and generate points exactly on the sphere.

    // We compute a point on a circle as follows:
    // Support you have a direction (x, y) and you want to find a
    // point (x', y') on a circle with rational radius r in that direction.
    // You do as follows:
    // 1) Compute phi = atan2(y,x) (using floating point).
    // 2) Compute tau = tan(phi/2) (using floating point).

    for (it = points.begin(); it != points.end(); ++it)
    {

      double phi = atan2(it->y(), it->x());
//      double tau = long(tan(phi/2) * 10) / double(10);
      FT sin, cos;
      get_exact_angles_on_circle(phi, cos, sin);
      out_points.push_back(Exact_point_2(10 * cos, 10 * sin));
    }
  }
  break;

  case 'o': // deterministic on circle with radius 10
  {
    Random_points_on_circle_2<Point, Creator> rnd_on_circ(1000);
    CGAL::copy_n(rnd_on_circ, n_points, std::back_inserter(points));

    for (int i = 0; i < n_points; ++i)
    {

      double phi = 360.0 * i / n_points;
//      double tau = long(tan(phi/2) * 10) / double(10);
      FT sin, cos;
      get_exact_angles_on_circle(phi, cos, sin);
      out_points.push_back(Exact_point_2(10 * cos, 10 * sin));
    }
  }
  break;

  case 'R': // random circles
  {
    for (std::size_t j = 0; j < n_points; ++j)
    {
      weights.push_back(
        CGAL::default_random.get_int(1,
                                     static_cast<int>(std::sqrt(n_points))
                                     * 100));
    }
    CGAL::copy_n(rnd_in_sq, n_points, std::back_inserter(points));
    copy_to_exact(points, out_points);
  }
    break;

  case 'A': // annulus
  {
    Random_points_on_circle_2<Point, Creator> rnd_on_circ(n_points);
    for (std::size_t j = 0; j < n_points; ++j)
    {
      weights.push_back(
        CGAL::default_random.get_int(n_points*5, n_points*20));
    }
    CGAL::copy_n(rnd_on_circ, n_points, std::back_inserter(points));
    CGAL::perturb_points_2 (points.begin(), points.end(), n_points*0.2);
    copy_to_exact(points, out_points);
  }
    break;

  case 'D': // degenerate annulus
  {
    Random_points_on_circle_2<Point, Creator> rnd_on_circ(n_points);
    FT outer_r = (int)(n_points);
    FT inner_r = (int)(n_points / 2);
    for (std::size_t j = 0; j < n_points; ++j)
    {
      weights.push_back(
        CGAL::default_random.get_int(n_points*5, n_points*20));

      Point p = *rnd_on_circ++;
      FT ecos, esin, edenom;
      CGAL::rational_rotation_approximation(FT(p.x()) / 10,
                                            FT(p.y()) / 10,
                                            esin, ecos, edenom,
                                            FT(1L),
                                            FT(int(10 * 100)));

      // get_exact_angles_on_circle(angle, ecos, esin, n_points * 10);

      FT scale = 0;
      if (j%2 == 0)
        scale = FT(outer_r) - FT(weights.back(), 100);
      else
        scale = FT(inner_r) + FT(weights.back(), 100);

      FT x = ecos * scale / edenom;
      FT y = esin * scale / edenom;
      out_points.push_back(Exact_point_2(x, y));
    }
  }
    break;


  case 'p': // random on parabola
    CGAL::copy_n(rnd_on_seg, n_points, std::back_inserter(points));
    for (Points_vector::iterator it=points.begin(); it!=points.end(); it++)
    {
      int slop = 1;
      *it = Point(it->x(), it->x() * it->x() * slop);
    }

    copy_to_exact(points, out_points);

  case 'S':
  {
    Random_points_on_circle_2<Point, Creator> rnd_on_circ(1000);
    for (std::size_t i = 0; i < n_points; ++i)
    {
      // we generate points on a sphere by first getting two angles
      // with rational coeff and then "rotating" the unit sphere.
      // we generate random angles using the generator of points on circle.
      Point p = *rnd_on_circ++;
      double phi = atan2(p.y(), p.x());
//      double tau = long(tan(phi/2) * 100) / double(100);
//      double tau = long(tan(phi/2) * 100) / double(100);

      FT sin_1, cos_1;
      get_exact_angles_on_circle(phi, sin_1, cos_1);

      p = *rnd_on_circ++;
      phi = atan2(p.y(), p.x());
//      tau = long(tan(phi/2) * 100) / double(100);

      FT sin_2, cos_2;
      get_exact_angles_on_circle(phi, sin_2, cos_2);

      // now that we have the two angles we generate the point.
      // rotate around the z-axis (1, 0, 0)
      FT x = cos_1;
      FT y = sin_1;
      FT z = 0;

      // roatate around x-axis
      FT old_y = y;
      y = old_y*cos_2;
      z = old_y*sin_2;

      out_points_on_sphere.push_back(Exact_point_3(x, y, z));
    }
  }
  break;

  case 'P':
  {
    boost::mt19937 boost_rand_gen;
    boost_rand_gen.seed(static_cast<unsigned int>(std::time(0)));
    boost::exponential_distribution<> distr(4);

    for (std::size_t i = 0; i < n_points; ++i)
    {
      double a = CGAL::default_random.get_double(-100, 100);
      double b = CGAL::default_random.get_double(-100, 100);
      double c = CGAL::default_random.get_double(-100, 100);

      // we want to have less circles near the origin.
      // we use exponential distribution to get more circles far from
      // the origin.
      boost::variate_generator<boost::mt19937&,
        boost::exponential_distribution<> > gen(boost_rand_gen, distr);

      // rat is in [-1, 1] with exponential distribution so that less
      // values are near the zero.
      double rat = 1 - gen() / 4;
      if (CGAL::default_random.get_bool() == true)
        rat = -rat;

      // The distance of a plane from the origin is:
      // $dist = |d| / \sqrt{a^2 + b^2 + c^2}$
      // so we set d as follows to get a plane which intersects the
      // unit sphere.
      double norm = sqrt(a*a + b*b + c*c);
      double d = rat * norm;

      // rounding
      FT A = int(a * 10);
      FT B = int(b * 10);
      FT C = int(c * 10);
      FT D = int(d * 10);

      out_planes_3.push_back(Exact_plane_3(A, B, C, D));
    }
  }
  break;

  case 'E':
  {

    // We wish to do the same as above, just that we want the planes to
    // intersect at (0, i_value, 0) --- we want the planes to intersect at a
    // single point.
    // This means that if the plane equation is: ax + by + cz + d = 0
    // then: b = -d/i_value
    //
    // If we plug this into the fact that we want the planes to intersect
    // the unit sphere we get (using the distance from the origin as above):
    // d = rat * \sqrt{a^2 + b^2 + (-d/i_value)^2}
    // which eventually becomes:
    // d = rat*\sqrt{\frac{a^2 + c^2}{1 - rat^2 / i_value^2}}
    //
    //
    const FT i_value = FT(4, 3); // should be more then 1.

    boost::mt19937 boost_rand_gen;
    boost_rand_gen.seed(static_cast<unsigned int>(std::time(0)));
    boost::exponential_distribution<> distr(5);

    for (std::size_t i = 0; i < n_points; ++i)
    {
      double a = CGAL::default_random.get_double(-100, 100);
      double c = CGAL::default_random.get_double(-100, 100);

      // we want to have less circles near the origin.
      // we use exponential distribution to get more circles far from
      // the origin.
      boost::variate_generator<boost::mt19937&,
        boost::exponential_distribution<> > gen(boost_rand_gen, distr);

      // rat is in [-1, 1] with exponential distribution so that less
      // values are near the zero.
      double rat = 1 - gen() / 5;
      if (CGAL::default_random.get_bool() == true)
        rat = -rat;

      // d = rat*\sqrt{\frac{a^2 + c^2}{1 - rat^2 / i_value^2}}
      double d = rat * sqrt((a*a + c*c) /
                            (1 - rat*rat/to_double(i_value*i_value)));

      // b = -d/i_value
      FT A = int(a * 10);
      FT C = int(c * 10);
      FT D = int(d * 10);
      FT B = -D/i_value;

      out_planes_3.push_back(Exact_plane_3(A, B, C, D));
    }
  }
  break;
  }

  CGAL_assertion((out_points.size() == n_points) ||           \
                 (out_points_on_sphere.size() == n_points) || \
                 (out_planes_3.size() == n_points));

  if (vm.count("sort"))
  {
    std::sort(out_points.begin(), out_points.end());
    std::sort(out_points_on_sphere.begin(), out_points_on_sphere.end());
    std::sort(out_planes_3.begin(), out_planes_3.end(), Less_then_planes());
  }
  else if(vm.count("shuffle"))
  {
    std::random_shuffle(out_points.begin(), out_points.end(), default_random);
    std::random_shuffle(out_points_on_sphere.begin(),
                        out_points_on_sphere.end(), default_random);
    std::random_shuffle(out_planes_3.begin(), out_planes_3.end(),
                        default_random);
  }

  if (out_points.empty() == false)
  {
    std::cout << out_points.size() << std::endl;
    Exact_point_2_vector::iterator begin_it = out_points.begin();
    for (Exact_point_2_vector::iterator it = begin_it;
         it != out_points.end(); ++it)
    {
      std::cout << *it;
      if (weights.empty() == false)
      {
        std::vector<int>::iterator wit = weights.begin();
        std::advance(wit, std::distance(begin_it, it));
        std::cout << " " << CGAL::Gmpq(*wit, 100);
      }
      std::cout << std::endl;
    }
  }

  if (out_points_on_sphere.empty() == false)
  {
    std::cout << out_points_on_sphere.size() << std::endl;
    Exact_point_3_vector::iterator begin_it = out_points_on_sphere.begin();
    for (Exact_point_3_vector::iterator it = begin_it;
         it != out_points_on_sphere.end(); ++it)
    {
      std::cout << *it;
      std::cout << std::endl;
    }
  }

  if (out_planes_3.empty() == false)
  {
    std::cout << out_planes_3.size() << std::endl;
    Exact_plane_3_vector::iterator begin_it = out_planes_3.begin();
    for (Exact_plane_3_vector::iterator it = begin_it;
         it != out_planes_3.end(); ++it)
    {
      std::cout << *it;
      std::cout << std::endl;
    }
  }

  return 0;
}
