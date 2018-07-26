#ifndef CGAL_LEVEL_OF_DETAIL_INTERNAL_UTILS
#define CGAL_LEVEL_OF_DETAIL_INTERNAL_UTILS

namespace CGAL {
namespace Level_of_detail {
namespace internal {

template <typename T>
class Indexer
{
  std::map<T, std::size_t> indices;
public:

  std::size_t operator() (const T& t)
  {
    return indices.insert (std::make_pair (t, indices.size())).first->second;
  }

  void clear() { indices.clear(); }
};

  
template<class FT, class SegmentsRandomAccessRange>
class Regular_segments_info_estimator {

public:

  using Segments = SegmentsRandomAccessRange;
  using Values = std::vector<FT>;

  Regular_segments_info_estimator(const Segments &segments) : 
    m_segments(segments), 
    m_length_threshold(-FT(1)) { 

    compute_length_threshold();
  }

  inline FT get_length_threshold() const {
                
    CGAL_precondition(m_length_threshold > FT(0));
    return m_length_threshold;
  }

  bool is_too_long_segment(const size_t segment_index) const {
                
    CGAL_precondition(m_length_threshold > FT(0));
    if (m_segments[segment_index]->get().squared_length() < m_length_threshold * m_length_threshold) return false;
    return true;
  }

private:
  const Segments &m_segments;
  FT m_length_threshold;

  void compute_length_threshold() {
                
    Values segment_lengths;
    compute_segment_lengths(segment_lengths);

    const FT mean = compute_mean(segment_lengths);
    const FT stde = compute_standard_deviation(segment_lengths, mean);

    const FT estimated_length_threshold = estimate_length_threshold(mean, stde);
    m_length_threshold = estimated_length_threshold;
  }

  void compute_segment_lengths(Values &segment_lengths) const {
    CGAL_precondition(m_segments.size() > 0);

    segment_lengths.clear();
    segment_lengths.resize(m_segments.size());

    for (size_t i = 0; i < m_segments.size(); ++i)
      segment_lengths[i] = compute_segment_length(i);
  }

  FT compute_segment_length(const size_t segment_index) const {
    return static_cast<FT>(CGAL::sqrt(CGAL::to_double(m_segments[segment_index]->get().squared_length())));
  }

  FT compute_mean(const Values &values) const {
    CGAL_precondition(values.size() > 0);

    FT mean = FT(0);
    for (size_t i = 0; i < values.size(); ++i)
      mean += values[i];
    mean /= static_cast<FT>(values.size());

    return mean;
  }

  FT compute_standard_deviation(const Values &values, const FT mean) const {
    CGAL_precondition(values.size() > 0);

    FT sum = FT(0);
    for (size_t i = 0; i < values.size(); ++i)
      sum += (values[i] - mean) * (values[i] - mean);
    sum /= static_cast<FT>(values.size());

    return static_cast<FT>(CGAL::sqrt(CGAL::to_double(sum)));
  }

  FT estimate_length_threshold(const FT mean, const FT stde) const {
    return mean + stde;
  }
};
  
template <typename InputParameters, typename InputSegments>
typename InputParameters::FT max_orientation_local
(const InputParameters& parameters, const InputSegments& segments, std::size_t segment_index)
{
  typedef typename InputParameters::FT FT;
  const FT value = parameters.max_angle_in_degrees();
  CGAL_precondition(value > FT(0));

  const Regular_segments_info_estimator<FT, InputSegments> segments_info_estimator(segments);
  if (segments_info_estimator.is_too_long_segment(segment_index))
    return parameters.small_fixed_orientation();

  return value;
}

template <typename InputParameters, typename InputSegments>
typename InputParameters::FT max_difference_local
(const InputParameters& parameters, const InputSegments& segments, std::size_t segment_index)
{
  typedef typename InputParameters::FT FT;
  const FT value = parameters.max_difference_in_meters();
  CGAL_precondition(value > FT(0));

  const Regular_segments_info_estimator<FT, InputSegments> segments_info_estimator(segments);
  if (segments_info_estimator.is_too_long_segment(segment_index))
    return parameters.small_fixed_difference();
                
  return value;
}

template <typename Point_3>
typename Kernel_traits<Point_3>::Kernel::Point_2
point_2_from_point_3 (const Point_3& point_3)
{
  return typename Kernel_traits<Point_3>::Kernel::Point_2 (point_3.x(), point_3.y());
}
  
template <typename Elements, typename PointMap>
typename Kernel_traits<typename boost::property_traits<PointMap>::value_type>::Kernel::Point_2
barycenter (const Elements& elements, PointMap point_map)
{
  typedef typename Kernel_traits<typename boost::property_traits<PointMap>::value_type>::Kernel Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  
  CGAL_precondition(elements.size() > 0);
                
  FT x = FT(0), y = FT(0), size = FT(0);

  for (typename Elements::const_iterator ce_it = elements.begin();
       ce_it != elements.end(); ++ce_it, size += FT(1))
  {
    const Point_2 &point = get (point_map, **ce_it);

    x += point.x();
    y += point.y();
  }

  x /= size;
  y /= size;

  return Point_2(x, y);
}

template <typename Kernel, typename FaceHandle>
typename Kernel::Point_2 barycenter (const FaceHandle& face_handle)
{
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  
  FT x = FT(0), y = FT(0);
  for (std::size_t i = 0; i < 3; ++i)
  {
    x += face_handle->vertex(i)->point().x();
    y += face_handle->vertex(i)->point().y();
  }

  x /= FT(3);
  y /= FT(3);

  return Point_2(x, y);
}

template<class Segment_2, class BoundingBox>
void compute_bounding_box_2(const std::vector<Segment_2>& segments,
                            BoundingBox& bounding_box) 
{
  typedef typename Kernel_traits<Segment_2>::Kernel Kernel;
  using FT        = typename Kernel::FT;
  using Point_2   = typename Kernel::Point_2;

  CGAL_precondition(segments.size() > 0);

  FT minx =  std::numeric_limits<FT>::max(), miny =  std::numeric_limits<FT>::max();
  FT maxx = -std::numeric_limits<FT>::max(), maxy = -std::numeric_limits<FT>::max();

  for (std::size_t i = 0; i < segments.size(); ++ i) {
    const Segment_2 &segment = segments[i];
                    
    const Point_2 &source = segment.source();
    const Point_2 &target = segment.target();

    minx = CGAL::min(minx, source.x()); minx = CGAL::min(minx, target.x());
    miny = CGAL::min(miny, source.y()); miny = CGAL::min(miny, target.y());

    maxx = CGAL::max(maxx, source.x()); maxx = CGAL::max(maxx, target.x());
    maxy = CGAL::max(maxy, source.y()); maxy = CGAL::max(maxy, target.y());
  }

  bounding_box.clear();
  bounding_box.push_back(Point_2(minx, miny));
  bounding_box.push_back(Point_2(maxx, miny));
  bounding_box.push_back(Point_2(maxx, maxy));
  bounding_box.push_back(Point_2(minx, maxy));
}

template<class Elements, class PointMap, class Plane_3, class BoundingBox>
void compute_bounding_box_3(const Elements& elements, const PointMap& point_map,
                            const Plane_3& plane, BoundingBox& bounding_box) 
{
  typedef typename Kernel_traits<typename boost::property_traits<PointMap>::value_type>::Kernel Kernel;
  using FT        = typename Kernel::FT;
  using Point_3   = typename Kernel::Point_3;
                
  CGAL_precondition(elements.size() > 0);

  FT minx =  std::numeric_limits<FT>::max(), miny =  std::numeric_limits<FT>::max();
  FT maxx = -std::numeric_limits<FT>::max(), maxy = -std::numeric_limits<FT>::max();
                
  FT z = FT(0), size = FT(0);
  for (typename Elements::const_iterator ce_it = elements.begin();
       ce_it != elements.end(); ++ce_it, size += FT(1)) {
					
    const Point_3 &point    = get(point_map, *ce_it);
    const Point_3 projected = plane.projection(point);

    minx = CGAL::min(minx, projected.x());
    miny = CGAL::min(miny, projected.y());

    maxx = CGAL::max(maxx, projected.x());
    maxy = CGAL::max(maxy, projected.y());

    z += projected.z();
  }
  z /= size;

  bounding_box.clear();
  bounding_box.push_back(Point_3(minx, miny, z));
  bounding_box.push_back(Point_3(maxx, miny, z));
  bounding_box.push_back(Point_3(maxx, maxy, z));
  bounding_box.push_back(Point_3(minx, maxy, z));
}

template <typename Plane_3, typename Point_2>
typename Kernel_traits<Point_2>::Kernel::Point_3
position_on_plane (const Plane_3& plane,
                   const Point_2& point)
{
  typedef typename Kernel_traits<Point_2>::Kernel Kernel;
  
  static typename Kernel::Vector_3 vertical (0., 0., 1.);
  typename Kernel::Line_3 line (typename Kernel::Point_3 (point.x(), point.y(), 0.), vertical);
    
  typename CGAL::cpp11::result_of<typename Kernel::Intersect_3(typename Kernel::Line_3, typename Kernel::Plane_3)>::type
    inter = CGAL::intersection (line, plane);
    
  if (inter)
    if (const typename Kernel::Point_3* p = boost::get<typename Kernel::Point_3>(&*inter))
      return *p;

  std::cerr << "Error: can't compute 3D position" << std::endl;
  return typename Kernel::Point_3 (0., 0., 0.);
}

template <typename Kernel>
struct Point_3_from_point_2_and_plane
{
  typedef typename Kernel::Point_2 argument_type;
  typedef typename Kernel::Point_3 result_type;

  const typename Kernel::Plane_3& plane;

  Point_3_from_point_2_and_plane (const typename Kernel::Plane_3& plane) : plane (plane) { }

  result_type operator() (const argument_type& a) const
  {
    return position_on_plane (plane, a);
  }
};
  
template <typename Kernel>
struct Segment_3_from_segment_2_and_plane
{
  typedef typename Kernel::Segment_2 argument_type;
  typedef typename Kernel::Segment_3 result_type;

  const typename Kernel::Plane_3& plane;

  Segment_3_from_segment_2_and_plane (const typename Kernel::Plane_3& plane) : plane (plane) { }

  result_type operator() (const argument_type& a) const
  {
    return result_type (position_on_plane (plane, a.source()),
                        position_on_plane (plane, a.target()));
  }
};
  
template <typename Triangulation>
void segment_semantic_faces (const Triangulation& triangulation,
                             std::vector<typename Triangulation::Face_handle>& ground_faces,
                             std::vector<typename Triangulation::Face_handle>& roof_faces,
                             std::vector<typename Triangulation::Face_handle>& vegetation_faces)
{
  for (typename Triangulation::Finite_faces_iterator
         it = triangulation.finite_faces_begin(); it != triangulation.finite_faces_end(); ++ it)
    if (it->info().visibility_label() == Visibility_label::INSIDE)
      roof_faces.push_back (it);
    else if (it->info().visibility_label() == Visibility_label::OUTSIDE)
      ground_faces.push_back (it);
}

template <typename Point_3, typename FaceHandle>
Point_3 point_3 (FaceHandle fh, std::size_t j)
{
  return Point_3 (fh->vertex(j)->point().x(), fh->vertex(j)->point().y(), fh->info().height(j));
}

template <typename Triangle_3, typename FaceHandle>
Triangle_3 triangle_3 (FaceHandle fh)
{
  typedef typename Kernel_traits<Triangle_3>::Kernel Kernel;
  typedef typename Kernel::Point_3 Point_3;
  return Triangle_3 (point_3<Point_3> (fh, 0), point_3<Point_3> (fh, 1), point_3<Point_3> (fh, 2));
}

} } } // namespace CGAL::Level_of_detail::internal

#endif
