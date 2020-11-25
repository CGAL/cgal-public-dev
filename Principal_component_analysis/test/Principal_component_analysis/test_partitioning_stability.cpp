#include <queue>

#include <CGAL/Bbox_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;
typedef Kernel::Line_3 Line_3;

static CGAL::Random rnd;

double distance (const Point_3& p, const Line_3& l)
{ return std::sqrt (CGAL::squared_distance (p, l)); }

template <typename Fitted>
void assert_quality (const std::vector<Point_3>& points, const Fitted& fitted)
{
  double mean_dist = 0;
  for (std::size_t i = 0; i < points.size(); ++ i)
  {
    double dist = distance (points[i], fitted);
    mean_dist += dist;
  }
  mean_dist /= points.size();

  std::cerr << "mean distance = " << mean_dist << std::endl;

  CGAL_assertion_code
    (double limit = 1e-4 * std::sqrt (CGAL::squared_distance (points.front(), points.back())));
  CGAL_assertion (mean_dist < limit);
}

void generate_random_points_on_line_around_bbox (
        const Line_3 &line,
        const CGAL::Bbox_3 &bbox,
        std::vector<Point_3>& target) {
  for(int i = 0; i < 100; i++) {
    double x = rnd.get_double(bbox.xmin(), bbox.xmax());
    double y = rnd.get_double(bbox.ymin(), bbox.ymax());
    double z = rnd.get_double(bbox.zmin(), bbox.zmax());
    Point_3 p = line.projection(Point_3(x,y,z));
    target.push_back(p);
  }
}

void test_partitioning_stability_segment()
{
  Point_3 p0 (rnd.get_double(), rnd.get_double(), rnd.get_double());
  Point_3 p1 (rnd.get_double(), rnd.get_double(), rnd.get_double());
  Segment_3 seg0(p0,p1);
  std::vector<Segment_3> initial_segment = {seg0};

  // Generating segment partitioning
  std::queue<Segment_3> queue;
  queue.push(seg0);
  int iteration = 0;
  while (iteration < 10) {
    Segment_3 seg = queue.front();
    queue.pop();
    CGAL::Random_points_on_segment_3<Point_3> g(seg[0],seg[1],rnd);
    std::vector<Point_3> rnd_pt;
    rnd_pt.reserve(1);
    std::copy_n( g, 1, std::back_inserter(rnd_pt));
    Segment_3 new_seg1(seg[0], rnd_pt[0]);
    Segment_3 new_seg2(rnd_pt[0], seg[1]);
    queue.push(new_seg1);
    queue.push(new_seg2);
    iteration++;
  }
  std::vector<Segment_3> partitioning;
  while (!queue.empty()) {
    partitioning.push_back(queue.front());
    queue.pop();
  }

  Line_3 initial_line;
  linear_least_squares_fitting_3(initial_segment.begin(),initial_segment.end(), initial_line, CGAL::Dimension_tag<1>());

  Line_3 partitioning_line;
  linear_least_squares_fitting_3(partitioning.begin(),partitioning.end(), partitioning_line, CGAL::Dimension_tag<1>());

  std::vector<Point_3> rnd_pt_on_partitioning_line;
  generate_random_points_on_line_around_bbox(partitioning_line,
                                             bbox_3(partitioning.begin(), partitioning.end()),
                                             rnd_pt_on_partitioning_line);
  assert_quality(rnd_pt_on_partitioning_line, initial_line);
}

int main()
{
  std::cerr << "Partitioning stability test with seed " << rnd.get_seed() << std::endl;

  test_partitioning_stability_segment();

  return EXIT_SUCCESS;
}
