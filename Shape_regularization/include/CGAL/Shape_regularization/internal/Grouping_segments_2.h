#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2

// #include <CGAL/license/Shape_regularization.h>

#include <vector>
#include <map>
#include <utility>

#include <CGAL/Shape_regularization/internal/Segment_data_2.h>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits,
    typename InputRange>
  class Grouping_segments_2 {
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using FT = typename GeomTraits::FT;
    using Segment_data = typename internal::Segment_data_2<Traits>;

    Grouping_segments_2(
      const InputRange& input_range) :
    m_input_range(input_range),
    m_theta_eps(FT(1) / FT(4)),
    m_tolerance(FT(1) / FT(1000000)) {

      CGAL_precondition(input_range.size() > 0);

    }

      void make_groups(const std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>> & t_ijs,
                       const std::map <std::pair<std::size_t, std::size_t>, std::pair<int, std::size_t>> & r_ijs,
                       const std::map <std::size_t, Segment_data> & segments,
                       const std::vector<FT> & orientations,
                       std::map<FT, std::vector<std::size_t>> & parallel_groups_by_angles) { 

      CGAL_precondition(orientations.size() > 0);
      parallel_groups_by_angles.clear();
      
      const std::size_t n = m_input_range.size();
      std::map <std::size_t, std::vector<std::size_t>> parallel_groups;
      std::map<std::size_t, int> segments_to_groups_hashmap;

      for (const auto & it : segments) {
        std::size_t seg_index = it.second.m_index;
        segments_to_groups_hashmap[seg_index] = -1;
      }

      std::size_t g = 0;
      int g_i, g_j;

      auto rel_it = r_ijs.begin();
      for (const auto & tar_it : t_ijs) {

        const std::size_t i = tar_it.first.first;
        const std::size_t j = tar_it.first.second;
        const std::size_t p = tar_it.second.second;
        CGAL_precondition(rel_it->second.second == p);
        const std::size_t r = rel_it->second.first;

        if (CGAL::abs(orientations[n + p]) < m_tolerance) { 

          if (segments_to_groups_hashmap[i] == -1 && segments_to_groups_hashmap[j] == -1) {
            switch (r) {
              case 0:
                // Then segments i and j belong to the same group of parallel segments.
                // We should create a group of segments, that is initialized with these two individuals.
                segments_to_groups_hashmap[i] = g;
                segments_to_groups_hashmap[j] = g;

                parallel_groups[g].push_back(i);
                parallel_groups[g].push_back(j);

                ++g;
                break;
              case 1:
                // The segments i and j are orthogonal.
                // We create two different groups of parallel segments.
                segments_to_groups_hashmap[i] = g;
                parallel_groups[g].push_back(i);
                segments_to_groups_hashmap[j] = ++g;
                parallel_groups[g].push_back(j);
                ++g;
                break;
            }
          }
          else if (segments_to_groups_hashmap[i] == -1 && segments_to_groups_hashmap[j] != -1) {
            switch (r) {
              case 0:
                // Then segment i is parallel to j, and can be assigned to the same group.
                g_j = segments_to_groups_hashmap[j];
                segments_to_groups_hashmap[i] = g_j;
                parallel_groups[g_j].push_back(i);
                break;
              case 1:
                // Then segment i is orthogonal to j, and we should initialize a new group with this segment.
                segments_to_groups_hashmap[i] = g;
                parallel_groups[g].push_back(i);
                ++g;
                break;
            }
          }
          else if (segments_to_groups_hashmap[i] != -1 && segments_to_groups_hashmap[j] == -1) {
            // Symmetrical situation to before.
            switch (r) {
              case 0:
                g_i = segments_to_groups_hashmap[i];
                segments_to_groups_hashmap[j] = g_i;
                parallel_groups[g_i].push_back(j);
                break;
              case 1:
                segments_to_groups_hashmap[j] = g;
                parallel_groups[g].push_back(j);
                ++g;
                break;
            }
          }
          else {
            g_i = segments_to_groups_hashmap[i];
            g_j = segments_to_groups_hashmap[j];
            if (g_i != g_j) {
              if (r == 0) {                       
                // Segments i and j have been assigned to different groups, but in fact
                // they are parallel and belong to the same group. That's why we merge them.
                for (const auto gr : parallel_groups[g_j]) {
                  segments_to_groups_hashmap[gr] = g_i;
                  parallel_groups[g_i].push_back(gr);
                }
                parallel_groups[g_j].clear();
              } else if (r == 1) {
                // We do nothing here.
              }
            }
          }
        }
        ++rel_it;
        
      }

      std::map<int, FT> angles;

      build_map_of_angles(orientations, segments, parallel_groups, 
                          segments_to_groups_hashmap, angles);

      // Try to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group.
      assign_segments_to_groups(segments, parallel_groups, segments_to_groups_hashmap, angles);

      build_parallel_groups_angle_map(segments_to_groups_hashmap, angles, parallel_groups_by_angles);

    } 

  private:
    const Input_range& m_input_range;
    std::vector<Segment_data> m_segments;
    const FT m_theta_eps;
    const FT m_tolerance;

  
    void build_map_of_angles(const std::vector<FT> & orientations,
                             const std::map <std::size_t, Segment_data> & segments,
                             std::map<std::size_t, std::vector<std::size_t>> & parallel_groups,
                             std::map<std::size_t, int> & segments_to_groups_hashmap,
                             std::map<int, FT> & angles) {

      for (const auto & sm_i : segments_to_groups_hashmap) {
        int g_i = sm_i.second;

        if (g_i != -1 && (angles.find(g_i) == angles.end())) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = segments.at(seg_index);
          FT theta = seg_data.m_orientation + orientations[seg_index];
          if (theta < FT(0)) 
            theta += FT(180);
          else if (theta > FT(180)) 
            theta -= FT(180);
          
          // Check if the angle that seems to be associated to this group of segments is not too close to another value.
          int g_j = -1;
          for (const auto & it_m : angles) {
            if (CGAL::abs(it_m.second - theta) < m_theta_eps) 
              g_j = it_m.first;
          }
          if (g_j == -1) 
            angles[g_i] = theta;
          else {                       
            // Merge groups.
            for (const auto gr : parallel_groups[g_i]) {
              segments_to_groups_hashmap[gr] = g_j;
              parallel_groups[g_j].push_back(gr);
            }
            parallel_groups[g_i].clear();
          }
        }
      }

    }

    void assign_segments_to_groups(const std::map <std::size_t, Segment_data> & segments,
                                   std::map<std::size_t, std::vector<std::size_t>> & parallel_groups,
                                   std::map<std::size_t, int> & segments_to_groups_hashmap,
                                   std::map<int, FT> & angles) {

      for (const auto & sm_i : segments_to_groups_hashmap) {
        int g_i = sm_i.second;
        if (g_i == -1) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = segments.at(seg_index);
          const FT alpha = seg_data.m_orientation;
          int g_j = -1;
          for (const auto & it_m : angles) {
            const FT alpha_j = it_m.second;
            for (int k = -1; k <= 1; ++k) {
              if (CGAL::abs(alpha_j - alpha + static_cast<FT>(k) * FT(180)) < m_theta_eps) {
                g_j = it_m.first;
                break;
              }
            }
            if (g_j != -1) 
              break;
          }
          if (g_j == -1) {   
            if (angles.size() > 0)               
              g_i = angles.rbegin()->first + 1;
            else
              g_i = 0;
            angles[g_i] = alpha;
          } 
          else 
            g_i = g_j;

          segments_to_groups_hashmap[seg_index] = g_i; // Segmentation fault (core dumped) for the example with 3 segments and bounds = 25
          parallel_groups[g_i].push_back(seg_index);
        }
      }

    }

    void build_parallel_groups_angle_map(const std::map<std::size_t, int> & segments_to_groups_hashmap,
                                         std::map<int, FT> & angles, 
                                         std::map<FT, std::vector<std::size_t>> & parallel_groups_by_angles) {

      for (const auto & it_m : angles) {
        const FT angle = it_m.second;
        if (parallel_groups_by_angles.find(angle) == parallel_groups_by_angles.end()) 
          parallel_groups_by_angles[angle] = std::vector<std::size_t>();
      }

      for (const auto & sm_i : segments_to_groups_hashmap) {
        // If segment s_i is included in a group of parallel segments,
        // then it should be assigned to a leaf of the regularization tree.
        const FT angle = angles.at(sm_i.second);
        if (parallel_groups_by_angles.find(angle) != parallel_groups_by_angles.end()) 
          parallel_groups_by_angles[angle].push_back(sm_i.first);
      }

    } 

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2