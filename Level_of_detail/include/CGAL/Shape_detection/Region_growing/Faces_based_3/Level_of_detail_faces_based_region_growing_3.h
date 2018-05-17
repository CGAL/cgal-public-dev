#ifndef CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_GROWING_H
#define CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_GROWING_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputMesh, class InputFaces>
		class Level_of_detail_planar_region_growing {

		public:
			typedef KernelTraits Kernel;
			typedef InputMesh 	 Mesh;
			typedef InputFaces 	 Faces;

			typedef typename Kernel::FT FT;

			typedef typename Kernel::Point_3  Point_3;
			typedef typename Kernel::Vector_3 Vector_3;

			typedef typename Mesh::Facet_const_handle    Facet_handle;
			typedef typename Mesh::Halfedge_const_handle Halfedge_handle;
			typedef typename Mesh::Vertex_const_handle   Vertex_handle;
			
			using Building_facets = std::vector<Facet_handle>;

			using Building_region  = std::vector<Facet_handle>;
			using Building_regions = std::vector<Building_region>;
			using Regions 		   = std::vector<Building_regions>;

			using Facet_map = std::map<Facet_handle, bool>;

			Level_of_detail_planar_region_growing(const Faces &building_faces) : 
			m_building_faces(building_faces), m_number_of_regions(-1), m_eps(FT(1) / FT(1000000)) { }

			void find_regions(Regions &regions) {
				assert(!m_building_faces.empty());

				regions.clear();
				regions.resize(m_building_faces.size());

				grow_regions(regions);
			}

			int get_number_of_regions() const {
				
				assert(m_number_of_regions >= 0);
				return m_number_of_regions;
			}

		private:
			const Faces &m_building_faces;

			int m_number_of_regions;
			const FT m_eps;

			void grow_regions(Regions &regions) {
				
				const size_t number_of_buildings = m_building_faces.size();
				assert(regions.size() == number_of_buildings);

				for (size_t i = 0; i < number_of_buildings; ++i) grow_building_regions(i, regions[i]);
				compute_number_of_regions(regions);
			}

			void compute_number_of_regions(const Regions &regions) {
				m_number_of_regions = 0;
				for (size_t i = 0; i < regions.size(); ++i)
					m_number_of_regions += static_cast<int>(regions[i].size());
			}

			void grow_building_regions(const size_t building_index, Building_regions &building_regions) {
				
				assert(building_index >= 0 && building_index < m_building_faces.size());
				building_regions.clear();

				const Building_facets &faces = m_building_faces[building_index];
				const size_t num_faces = faces.size();

				std::map<Facet_handle, bool> handled_faces;
				for (size_t i = 0; i < num_faces; ++i) {
					
					const Facet_handle &fh = faces[i];
					handled_faces[fh] = false;
				}

				for (size_t i = 0; i < num_faces; ++i) {
					
					const Facet_handle &fh = faces[i];
					if (!handled_faces.at(fh)) {
						
						Building_region new_region;
						grow_new_region(fh, faces, new_region, handled_faces);
						building_regions.push_back(new_region);
					}
				}
			}

			void grow_new_region(const Facet_handle &fh, const Building_facets &faces, Building_region &building_region, Facet_map &handled_faces) {

				building_region.push_back(fh);
				handled_faces[fh] = true;

				const size_t num_faces = faces.size();
				for (size_t i = 0; i < num_faces; ++i) grow(fh, faces[i], building_region, handled_faces);
			}

			void grow(const Facet_handle &curr, const Facet_handle &next, Building_region &building_region, Facet_map &handled_faces) {

				if (handled_faces.at(next)) return;
				if (!should_be_added(curr, next)) return;

				building_region.push_back(next);
				handled_faces[next] = true;
			}

			bool should_be_added(const Facet_handle &curr, const Facet_handle &next) {

				Vector_3 n1 = get_normal(curr);
				Vector_3 n2 = get_normal(next);

				const FT scalar = CGAL::scalar_product(n1, n2);
				if (CGAL::abs(scalar - FT(1)) < m_eps || CGAL::abs(scalar + FT(1)) < m_eps)
					if (are_coplanar(curr, next))
						return true;

				return false;
			}

			Vector_3 get_normal(const Facet_handle &fh) {

				Halfedge_handle he = fh->halfedge();
				const Point_3 p1 = he->vertex()->point();

				he = he->next();
				const Point_3 p2 = he->vertex()->point();

				he = he->next();
				const Point_3 p3 = he->vertex()->point();

				const Vector_3 v1 = Vector_3(p1, p2);
				const Vector_3 v2 = Vector_3(p1, p3);

				const Vector_3 cross = CGAL::cross_product(v1, v2);
				const Vector_3 n = cross / static_cast<FT>(CGAL::sqrt(CGAL::to_double(cross.squared_length())));

				return n;
			}

			bool are_coplanar(const Facet_handle &curr, const Facet_handle &next) {

				std::vector<Point_3> v1;
				get_vertices(curr, v1);

				std::vector<Point_3> v2;
				get_vertices(next, v2);

				int count = 0;
				for (size_t j = 0; j < v2.size(); ++j)	
					if (is_coplanar(v1[0], v1[1], v1[2], v2[j])) ++count;

				if (count >= 1) return true;
				return false;
			}

			void get_vertices(const Facet_handle &fh, std::vector<Point_3> &v) {
				Halfedge_handle he = fh->halfedge();

				const Point_3 p1 = he->vertex()->point();

				he = he->next();
				const Point_3 p2 = he->vertex()->point();

				he = he->next();
				const Point_3 p3 = he->vertex()->point();

				v.push_back(p1);
				v.push_back(p2);
				v.push_back(p3);
			}

			bool is_coplanar(const Point_3 &p1, const Point_3 &p2, const Point_3 &p3, const Point_3 &p4) {

				const Vector_3 cross = CGAL::cross_product(Vector_3(p1, p2), Vector_3(p1, p4));
				const FT result = CGAL::scalar_product(cross, Vector_3(p1, p3));

				if (result < m_eps) return true;
				return false;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_GROWING_H