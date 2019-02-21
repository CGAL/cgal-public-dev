#ifndef CGAL_LEVELS_OF_DETAIL_STRUCTURES_H
#define CGAL_LEVELS_OF_DETAIL_STRUCTURES_H

// STL includes.
#include <vector>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  struct Polygon_face_2 {

  public:
    using Traits = GeomTraits;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    std::vector<Point_2> vertices;
    Visibility_label visibility;
    std::vector<std::size_t> neighbors;
    std::vector<Segment_2> edges;

    Polygon_face_2() :
    visibility(Visibility_label::OUTSIDE) 
    { }
  };

  template<typename GeomTraits>
  struct Polyhedron_facet_3 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    std::vector<Point_3> vertices;
    std::vector< std::vector<std::size_t> > faces;

    Visibility_label visibility;
    std::vector<std::size_t> neighbors;

    FT weight = FT(0);
    FT in = FT(0);
    FT out = FT(0);

    Polyhedron_facet_3() :
    visibility(Visibility_label::OUTSIDE) 
    { }
  };

  template<typename GeomTraits>
  struct Graphcut_face_3 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    using Data = std::pair<int, int>;
		using Pair = std::pair<Data, Data>;

		Data data;
		Pair neighbors;

		FT weight  = -FT(1);
		FT quality = -FT(1);

    bool is_valid = true;
  };

  template<typename GeomTraits>
  struct Face_info {

  public:
    using Traits = GeomTraits;
  };

  template<typename GeomTraits>
  struct Vertex_info {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    FT height = -FT(100000000000000);
  };

  template<typename GeomTraits>
  struct Roof {

  public:
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;

    Roof(std::vector<Point_3> roof) :
    vertices(roof)
    { }

    std::vector<Point_3> vertices;
  };

  template<typename GeomTraits>
  struct Wall {

  public:
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;

    Wall(std::vector<Point_3> wall) :
    vertices(wall)
    { }

    std::vector<Point_3> vertices;
  };

  template<typename GeomTraits>
  struct Building {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Triangle_2 = typename Traits::Triangle_2;
    using Vector_3 = typename Traits::Vector_3;
    using Point_3 = typename Traits::Point_3;
    
    FT height = FT(0);
    std::vector<Triangle_2> footprint;
    std::vector<Segment_2> boundaries;

    std::size_t cluster_index;
    std::vector< std::vector<std::size_t> > roof_indices;
    std::vector<Vector_3> normals;
    std::vector< std::vector<Point_3> > approximate_roofs;
    std::vector< std::vector<Point_3> > approximate_walls;
    std::vector<Point_3> approximate_ground;

    std::vector< Polyhedron_facet_3<Traits> > polyhedrons;
    std::vector< Graphcut_face_3<Traits> > graphcut_faces;

    std::vector< Roof<Traits> > roofs;
    std::vector< Wall<Traits> > walls;
  };

  template<typename GeomTraits>
  struct Tree {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Triangle_2 = typename Traits::Triangle_2;

    FT radius = FT(1);
    Point_2 center;

    FT height = FT(0);
    std::vector<Triangle_2> footprint;
    std::vector<Segment_2> boundaries;

    std::size_t cluster_index;
  };

  template<typename GeomTraits>
  struct Smooth_ground {

  public:
    using Traits = GeomTraits;
    using Triangle_3 = typename Traits::Triangle_3;

    std::vector<Triangle_3> triangles;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_STRUCTURES_H
