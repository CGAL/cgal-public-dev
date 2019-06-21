#ifndef CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2
#define CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2

// #include <CGAL/license/Shape_regularization.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits, 
    typename InputRange, 
    typename SegmentMap> //property map - read!
  class Delaunay_neighbor_query_2 {

  public:

    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using FT = typename GeomTraits::FT;
    using Point = typename GeomTraits::Point_2;
    using Segment = typename GeomTraits::Segment_2;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<int, GeomTraits>;
    using DS = CGAL::Triangulation_data_structure_2<VB>;
    using DT = CGAL::Delaunay_triangulation_2<GeomTraits, DS>;

    using Vertex_iterator = typename DT::Finite_vertices_iterator;
    using Vertex_circulator = typename DT::Vertex_circulator;


    Delaunay_neighbor_query_2(
      const InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) {

      CGAL_precondition(input_range.size() > 0);

      build_delaunay_triangulation();

    }

    void operator()(int i, std::vector<int> & result) { 

      // returns std::vector indicies of neighbors
      // Uses Delaunay triangulation to find neighbors

      for (Vertex_iterator vit = m_dt.finite_vertices_begin(); vit != m_dt.finite_vertices_end(); ++vit) {
        // std::cout << vit->point() << " -> info() = " << vit->info() << std::endl; // for debugging purposes
        if(vit->info() == i) {
          Vertex_circulator vc(vit);
          do {
            if(!m_dt.is_infinite(vc)) {
              result.push_back(vc->info());
            }
            --vc;
          } while (vc != m_dt.incident_vertices(vit));
      // incident verticies (Vertex_circulator 	incident_vertices (Vertex_handle v) const)
      // takes vit and returns vetrex_circulator
      // use do {} while() 
      // use  () to check if the vertex is connected to the original one
      // if compiler won't like Vertex_iterator use "static_cast" to convert types
          return;
        }
      }

    }

    //should be very similar to the K_neighbor_query class
    // in constructor:
    // 1) Check all segments in a "for" loop
    // 2) Call SegmentMap which returns a segment (center)
    // 3) Calculate a sentral point for each segment
    // 4) Send central points from each segment to Delaunay triangulation 
    // (from the old code); read Delaunay triangulation documentation
    // 5) Save Delaunay triangulation

// find neighbours in operator: give input index, return one ring heignhborhood of neighbours.
  private:

    // Fields.
    const Input_range& m_input_range;
    const Segment_map  m_segment_map;
    DT                 m_dt;

    Point compute_barycentre(const Segment& m_segment) {
      const FT half = FT(1) / FT(2);
      const Point &source = m_segment.source();
      const Point &target = m_segment.target();
      const FT x = half * (source.x() + target.x());
      const FT y = half * (source.y() + target.y());
      return Point(x, y);
    }

    void build_delaunay_triangulation() {
      m_dt.clear();
      for(int i = 0; i < m_input_range.size(); ++i) {
        m_dt.insert(compute_barycentre(m_input_range[i]))->info() = i;
      }
    }


  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2