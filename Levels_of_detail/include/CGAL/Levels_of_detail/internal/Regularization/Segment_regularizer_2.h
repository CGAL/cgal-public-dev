#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <memory>
#include <utility>

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/internal/Fitting/Segment_to_points_fitter.h>

#include <CGAL/Level_of_detail/internal/Regularization/Segment_regularizer_parameters.h>

#include <CGAL/Level_of_detail/internal/Regularization/Angles_regularizer.h>
#include <CGAL/Level_of_detail/internal/Regularization/Ordinates_regularizer.h>

#include <CGAL/Level_of_detail/internal/Regularization/Regular_segment.h>

namespace CGAL {

namespace Level_of_detail {


template<class GeomTraits>
class Segment_regularizer_2 {

public:
  using Kernel = GeomTraits;

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Segment_2 = typename Kernel::Segment_2;

  using Segment_fitter      = Segment_to_points_fitter<Kernel>;
  
  using Regular_segment     = Regular_segment<Kernel>;
  using Regular_segments    = std::vector<Regular_segment *>;
  using In_regular_segments = std::vector<Regular_segment>;
            
  using Angles_regularizer    = Angles_regularizer<Kernel>;
  using Ordinates_regularizer = Ordinates_regularizer<Kernel>;

  using Parameters = Segment_regularizer_parameters<FT>;
  using Tree       = typename Angles_regularizer::Tree;

  Segment_regularizer_2(const Parameters &parameters) : 
    m_parameters(parameters)
  { }

  ~Segment_regularizer_2() {
    m_input_segments.clear();
  }

  
  void regularize(const std::vector<std::vector<std::size_t> >& elements,
                  const std::vector<Point_2>& points,
                  std::vector<Segment_2>& output) {
    if (elements.size() == 0) return;

    // Copy input segments in the internal storage.
    create_input(elements, points);
                
    // Make segments parallel and orthogonal.
    if (m_parameters.optimize_angles()) 
      apply_angles_regularization();

    // Make segments collinear.
    if (m_parameters.optimize_angles() && m_parameters.optimize_ordinates()) 
      apply_ordinates_regularization();

    // Return all resulting segments.
    output.reserve (m_input_segments.size());
    for (size_t i = 0; i < m_input_segments.size(); ++i)
      output.push_back (m_input_segments[i]->get());
  }

private:            
  const Parameters &m_parameters;
            
  Regular_segments    m_input_segments;
  In_regular_segments m_in_segments;

  std::shared_ptr<Angles_regularizer>    m_angles_regularizer_ptr;
  std::shared_ptr<Ordinates_regularizer> m_ordinates_regularizer_ptr;
            
  void create_input(const std::vector<std::vector<std::size_t> >& elements,
                    const std::vector<Point_2>& points) {
                
    CGAL_precondition(elements.size() > 0);

    m_in_segments.clear();
    m_in_segments.resize(elements.size());

    for (std::size_t i = 0; i < elements.size(); ++ i) {
      Segment_2 segment;
      Segment_fitter().fit_segment_2 (elements[i], points, segment);
      m_in_segments[i] = Regular_segment(i, segment);
    }

    m_input_segments.clear();
    m_input_segments.resize(m_in_segments.size());

    for (size_t i = 0; i < m_in_segments.size(); ++i)
      m_input_segments[i] = &m_in_segments[i];
  }

  void apply_angles_regularization() {

    m_angles_regularizer_ptr = std::make_shared<Angles_regularizer>(m_input_segments, m_parameters);
    m_angles_regularizer_ptr->regularize();
  }

  void apply_ordinates_regularization() {
    Tree *tree_pointer = m_angles_regularizer_ptr->get_tree_pointer();

    m_ordinates_regularizer_ptr = std::make_shared<Ordinates_regularizer>(m_input_segments, tree_pointer, m_parameters);
    m_ordinates_regularizer_ptr->regularize();
  }

  template<class Segment_map, class Output>
  void create_output(const Segment_map &segment_map, Output &output) {

    for (size_t i = 0; i < m_input_segments.size(); ++i)
      put(segment_map, output, m_input_segments[i]->get());
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H
