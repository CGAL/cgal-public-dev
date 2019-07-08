#ifndef CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2
#define CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2

// #include <CGAL/license/Shape_regularization.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>
#include <CGAL/Shape_regularization/internal/utils.h>

#include <vector>


namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits, 
    typename InputRange, 
    typename SegmentMap>
  class Delaunay_neighbor_query_2 {

  public:

    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using FT = typename GeomTraits::FT;
    using Point = typename GeomTraits::Point_2;
    using Segment = typename GeomTraits::Segment_2;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, GeomTraits>;
    using DS = CGAL::Triangulation_data_structure_2<VB>;
    using DT = CGAL::Delaunay_triangulation_2<GeomTraits, DS>;

    using Vertex_circulator = typename DT::Vertex_circulator;


    /* n constructor:
    1) Check all segments in a "for" loop
    2) Call SegmentMap which returns a segment (center)
    3) Calculate a sentral point for each segment
    4) Send central points from each segment to Delaunay triangulation 
    (from the old code); read Delaunay triangulation documentation
    5) Save Delaunay triangulation */

    Delaunay_neighbor_query_2(
      InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) {

      CGAL_precondition(input_range.size() > 0);

      build_delaunay_triangulation();
      build_map_of_neighbours();

    }

    /* returns std::vector indicies of neighbors
    Uses Delaunay triangulation to find neighbors */
    void operator()(std::size_t i, std::vector<std::size_t> & neighbors) { 

      neighbors.clear();
      neighbors = m_map_of_neighbours[i];

    }

  private:

    // Fields.
    Input_range& m_input_range;
    const Segment_map  m_segment_map;
    DT                 m_dt;
    std::vector <std::vector<std::size_t>> m_map_of_neighbours;

    void build_delaunay_triangulation() {
      m_dt.clear();
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const Segment& seg = get(m_segment_map, *(m_input_range.begin() + i));
        const Point& source = seg.source();
        const Point& target = seg.target();
        const Point middle_point = internal::compute_middle_point(source, target);
        auto vh = m_dt.insert(middle_point);
        vh->info() = i;
      }
    }

    void build_map_of_neighbours() {
      m_map_of_neighbours.clear();
      m_map_of_neighbours.resize(m_input_range.size());
      for (auto vit = m_dt.finite_vertices_begin(); vit != m_dt.finite_vertices_end(); ++vit) {
        Vertex_circulator vc(vit);
        do {
          if(!m_dt.is_infinite(vc)) {
            CGAL_precondition(vit->info() >= 0 && vit->info() < m_input_range.size() 
              && vc->info() >= 0 && vc->info() < m_input_range.size());
            m_map_of_neighbours[vit->info()].push_back(vc->info());
          }
          --vc;
        } while (vc != m_dt.incident_vertices(vit));
      } 

    }


  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2
