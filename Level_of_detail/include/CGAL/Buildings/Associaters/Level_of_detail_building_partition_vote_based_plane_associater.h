#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_VOTE_BASED_PLANE_ASSOCIATER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_VOTE_BASED_PLANE_ASSOCIATER_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Associaters/Level_of_detail_building_partition_naive_plane_associater.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputBuilding>
		class Level_of_detail_building_partition_vote_based_plane_associater {
            
        public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef InputBuilding  Building;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_3 = typename Kernel::Triangle_3;

            using Roof              = typename Building::Roof;
            using Envelope_input    = typename Building::Data_triangles;
            using Associated_planes = typename Roof::Associated_planes;
            
            using Boundary = std::vector<Point_3>;
            using Log      = CGAL::LOD::Mylog;

            using Naive_plane_associater = CGAL::LOD::Level_of_detail_building_partition_naive_plane_associater<Kernel, Building>;

            using Shape_votes = std::vector<FT>;
            using Roof_votes  = std::vector<Shape_votes>;

            using Index   = int;
			using Indices = std::vector<Index>;

            using Polygon  = std::vector<Point_2>; 
            using Polygons = std::vector<Polygon>;

            using Coordinates = std::vector<FT>;

            using Mean_value = CGAL::Barycentric_coordinates::Mean_value_2<Kernel>;
            using Mean_value_coordinates = CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel>;

            Level_of_detail_building_partition_vote_based_plane_associater(const Input &input, const Building &building, const FT reference_height) :
            m_input(input),
            m_building(building),
            m_reference_height(reference_height) { 

                set_initial_roof_votes();
                create_polygons();
                compute_final_roof_votes();
            }

            void find_associated_planes(const size_t roof_face_index, bool &is_plane_index, Associated_planes &associated_planes) const {
                if (naive_should_be_used(roof_face_index, associated_planes)) return;

                associated_planes.clear();
                assert(m_roof_votes[roof_face_index].size() != 0);
                
                FT maxv = m_roof_votes[roof_face_index][0];
                size_t plane_index = 0;

                for (size_t i = 1; i < m_roof_votes[roof_face_index].size(); ++i) {
                    const FT value = m_roof_votes[roof_face_index][i];

                    if (value > maxv) {

                        maxv = value;
                        plane_index = i;
                    }
                }
                    
                is_plane_index = true;
                associated_planes.push_back(plane_index);
            }

        private:
            const Input    &m_input;
            const Building &m_building;
            
            const FT m_reference_height;
            Polygons m_polygons;

            Roof_votes             m_roof_votes;
            Naive_plane_associater m_naive_plane_associater;

            void set_initial_roof_votes() {

                const size_t num_roof_faces  = m_building.roofs.size();
                const size_t num_roof_shapes = m_building.shapes.size();

                assert(num_roof_faces  != 0);
                assert(num_roof_shapes != 0);

                m_roof_votes.clear();
                m_roof_votes.resize(num_roof_faces);

                for (size_t i = 0; i < num_roof_faces; ++i) {
                    
                    m_roof_votes[i].clear();
                    m_roof_votes[i].resize(num_roof_shapes, FT(0));
                }
            }

            void create_polygons(){

                const size_t num_roof_faces = m_building.roofs.size();
                assert(num_roof_faces != 0);

                m_polygons.clear();
                m_polygons.resize(num_roof_faces);

                for (size_t i = 0; i < num_roof_faces; ++i) {
                    const size_t num_points = m_building.roofs[i].boundary.size();

                    m_polygons[i].clear();
                    m_polygons[i].resize(num_points);

                    for (size_t j = 0; j < num_points; ++j) {
                        
                        const Point_3 &p = m_building.roofs[i].boundary[j];
                        m_polygons[i][j] = Point_2(p.x(), p.y());
                    }
                }
            }

            // This function can be significantly accelerated taking into account the proximity of points inside the given shape!
            void compute_final_roof_votes() {
                const auto &shapes = m_building.shapes;

                const size_t num_roof_shapes = shapes.size();
                assert(num_roof_shapes != 0);

                for (size_t i = 0; i < num_roof_shapes; ++i) {
                    const Indices &indices = shapes[i];

                    for (size_t j = 0; j < indices.size(); ++j) {
                        
                        const Point_3 &p = m_input.point(indices[j]);
                        const int roof_face_index = find_corresponding_roof_face(p);

                        if (roof_face_index >= 0 && roof_face_index < m_roof_votes.size())
                            m_roof_votes[roof_face_index][i] += FT(1);
                    }
                }
            }

            int find_corresponding_roof_face(const Point_3 &query) const {
                const Point_2 p = Point_2(query.x(), query.y());

                for (size_t i = 0; i < m_polygons.size(); ++i) {
                    Mean_value_coordinates mean_value_coordinates(m_polygons[i].begin(), m_polygons[i].end());

                    Coordinates coordinates;
                    mean_value_coordinates(p, std::back_inserter(coordinates));

                    bool found = true;
                    for (size_t j = 0; j < coordinates.size(); ++j) {
                        if (coordinates[j] < FT(0)) {
                            
                            found = false;
                            break;
                        }
                    }
                    if (found) return i;
                }
                return -1;
            }

            bool naive_should_be_used(const size_t roof_face_index, Associated_planes &associated_planes) const {

                for (size_t i = 0; i < m_roof_votes[roof_face_index].size(); ++i)
                    if (m_roof_votes[roof_face_index][i] != FT(0)) return false;

                m_naive_plane_associater.find_associated_planes(m_building.roofs[roof_face_index].boundary, m_building.envelope_input, m_reference_height, associated_planes);
                return true;
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_VOTE_BASED_PLANE_ASSOCIATER_H