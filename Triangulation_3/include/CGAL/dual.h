#include <utility>
#include <iterator>
#include <vector>
#include <map>

#include <boost/bind.hpp>
#include <boost/function_output_iterator.hpp>

template < class Delaunay_triangulation_3, class Polyhedron, class Polyhedron_incremental_builder, class Bbox >
Polyhedron
dual(const Delaunay_triangulation_3& t, typename Delaunay_triangulation_3::Vertex_handle v, const Bbox& bbox/* = Bbox(0., 0., 0., 0., 0., 0.)*/)
{
  //TODO: sweep typedefs
  typedef typename Delaunay_triangulation_3::Point Point;
  typedef typename Delaunay_triangulation_3::Vertex Vertex;
  typedef typename Delaunay_triangulation_3::Edge Edge;
  typedef typename Delaunay_triangulation_3::Facet Facet;
  typedef typename Delaunay_triangulation_3::Cell Cell;
  typedef typename Delaunay_triangulation_3::Vertex_handle Vertex_handle;
  typedef typename Delaunay_triangulation_3::Cell_handle Cell_handle;
  typedef typename Delaunay_triangulation_3::Cell_circulator Cell_circulator;

  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_precondition(!t.is_infinite(v) );
  CGAL_triangulation_precondition(t.dimension() == 3);
  CGAL_triangulation_expensive_precondition(t.tds().is_vertex(v) );

  typedef std::map<Cell_handle, size_t> Cell_id_map;
  // Indices correspond to the incident to v cells and to its dual vertices.
  Cell_id_map cell_ids;
  {
    typedef std::pair<typename Cell_id_map::iterator, bool> Iter_bool_pair;
    typedef typename Cell_id_map::value_type Value_type;
    typedef Iter_bool_pair(Cell_id_map::*Insert)(const Value_type&);
    t.incident_cells(v, boost::make_function_output_iterator(
      boost::bind(static_cast<Insert>(&Cell_id_map::insert), &cell_ids,
        // Need const references for C++11 compability.
        boost::bind(&std::make_pair<const Cell_handle&, const size_t&>,
          _1, boost::bind(&Cell_id_map::size, &cell_ids)))));
  }

  // Dual vertices of finite cells incident to v.
  std::vector<Point> points(cell_ids.size());
  // A first found facet which is incident to v and which belongs to the hull.
  // of the triangulation, if any.
  Facet first_hull_facet;
  bool hull_facet_found = false;
  for(typename Cell_id_map::iterator cit = cell_ids.begin(), end = cell_ids.end(); cit != end; ++cit) {
    Cell_handle cell = cit->first;
    int inf_v_id;
    if(cell->has_vertex(t.infinite_vertex(), inf_v_id) && !hull_facet_found) {
      hull_facet_found = true;
      first_hull_facet = t.mirror_facet(std::make_pair(cell, inf_v_id));
    }
    points[cit->second] = t.dual(cell);
  }

  // First step: add finite facets to the resulting polyhedron.
  std::vector<Edge> edges;
  t.incident_edges(v, std::back_inserter(edges));

  Polyhedron result;
  Polyhedron_incremental_builder builder(result.hds()/*TODO: TMP*/, true);
  builder.begin_surface(cell_ids.size(), edges.size(), 3 * cell_ids.size());
  for(typename std::vector<Point>::iterator pit = points.begin(), pend = points.end(); pit != pend; ++pit) {
    builder.add_vertex(*pit);
  }

  for(typename std::vector<Edge>::iterator eit = edges.begin(), eend = edges.end(); eit != eend; ++eit) {
    Cell_circulator ccir = t.incident_cells(*eit);
    Cell_circulator cend = ccir;
    builder.begin_facet();
    do {
      builder.add_vertex_to_facet(cell_ids[ccir]);
      ++ccir;
    } while ( ccir != cend );
    builder.end_facet();
  }
  builder.end_surface();

  if(!hull_facet_found) {
  // Dual cell is finite, it does not need to be intersected with the bbox.
    CGAL_triangulation_expensive_postcondition(result.is_valid());
    CGAL_triangulation_postcondition(result.is_closed());

    return result;
  }

  // Second step: close the resulting infinite polyhedron, if any.
  // TODO: implement
  return result;
}