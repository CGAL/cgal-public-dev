#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_ORDINATES_2
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_ORDINATES_2

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
  class Grouping_segments_ordinates_2 {
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using FT = typename GeomTraits::FT;
    using Segment_data = typename internal::Segment_data_2<Traits>;

    Grouping_segments_ordinates_2(
      const InputRange& input_range) :
    m_input_range(input_range),
    m_y_eps(FT(1)),
    m_tolerance(FT(1) / FT(1000000)) {

      CGAL_precondition(input_range.size() > 0);

    }

    void make_groups(const std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>> & t_ijs,
                     const std::map <std::size_t, Segment_data> & segments,
                     const std::vector<FT> & qp_result,
                     std::map<FT, std::vector<std::size_t>> & groups_by_value,
                     const std::map <std::pair<std::size_t, std::size_t>, std::pair<int, std::size_t>> & r_ijs = 
                     std::map <std::pair<std::size_t, std::size_t>, std::pair<int, std::size_t>>()) { 
      
      CGAL_precondition(qp_result.size() > 0);
      groups_by_value.clear();
      
      const std::size_t n = m_input_range.size();
      std::map<std::size_t, std::vector<std::size_t>> groups;
      std::map<std::size_t, int> segments_to_groups_hashmap;

      for (const auto & it : segments) {
        std::size_t seg_index = it.second.m_index;
        segments_to_groups_hashmap[seg_index] = -1;
      }

      build_initial_groups(n, t_ijs, r_ijs, segments, qp_result, groups, segments_to_groups_hashmap);

      std::map<int, FT> values;

      build_map_of_values(qp_result, segments, groups, 
                          segments_to_groups_hashmap, values);

      // Try to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group.
      assign_segments_to_groups(segments, groups, segments_to_groups_hashmap, values);

      build_groups_by_value(segments_to_groups_hashmap, values, groups_by_value);

    }

  private:
    const Input_range& m_input_range;
    const FT m_y_eps;
    const FT m_tolerance;

    void build_initial_groups(const std::size_t n,
                              const std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>> & t_ijs,
                              const std::map <std::pair<std::size_t, std::size_t>, std::pair<int, std::size_t>> & r_ijs,
                              const std::map <std::size_t, Segment_data> & segments,
                              const std::vector<FT> & qp_result,
                              std::map<std::size_t, std::vector<std::size_t>> & groups,
                              std::map<std::size_t, int> & segments_to_groups_hashmap) {

      std::size_t g = 0;
      int g_i = -1; 
      int g_j = -1;

      auto rel_it = r_ijs.begin();
      for (const auto & tar_it : t_ijs) {

        const std::size_t i = tar_it.first.first;
        const std::size_t j = tar_it.first.second;
        const std::size_t p = tar_it.second.second;
        int r = 0;
        if (rel_it != r_ijs.end()) {
          CGAL_precondition(rel_it->second.second == p);
          r = rel_it->second.first;
        }

        if (CGAL::abs(qp_result[n + p]) < m_tolerance) { 

          if (segments_to_groups_hashmap[i] == -1 && segments_to_groups_hashmap[j] == -1) {
            switch (r) {
              case 0:
                // Then segments i and j belong to the same group of parallel segments.
                // We should create a group of segments, that is initialized with these two individuals.
                segments_to_groups_hashmap[i] = g;
                segments_to_groups_hashmap[j] = g;

                groups[g].push_back(i);
                groups[g].push_back(j);

                ++g;
                break;
              case 1:
                // The segments i and j are orthogonal.
                // We create two different groups of parallel segments.
                segments_to_groups_hashmap[i] = g;
                groups[g].push_back(i);
                segments_to_groups_hashmap[j] = ++g;
                groups[g].push_back(j);
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
                groups[g_j].push_back(i);
                break;
              case 1:
                // Then segment i is orthogonal to j, and we should initialize a new group with this segment.
                segments_to_groups_hashmap[i] = g;
                groups[g].push_back(i);
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
                groups[g_i].push_back(j);
                break;
              case 1:
                segments_to_groups_hashmap[j] = g;
                groups[g].push_back(j);
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
                for (const auto gr : groups[g_j]) {
                  segments_to_groups_hashmap[gr] = g_i;
                  groups[g_i].push_back(gr);
                }
                groups[g_j].clear();
              } else if (r == 1) {
                // We do nothing here.
              }
            }
          }
        }
        if(rel_it != r_ijs.end())
          ++rel_it;
      }

    }

    void build_map_of_values(const std::vector<FT> & qp_result,
                             const std::map <std::size_t, Segment_data> & segments,
                             std::map<std::size_t, std::vector<std::size_t>> & groups,
                             std::map<std::size_t, int> & segments_to_groups_hashmap,
                             std::map<int, FT> & values) {

      for (const auto & sm_i : segments_to_groups_hashmap) {
        const int g_i = sm_i.second;

        if (g_i != -1 && (values.find(g_i) == values.end())) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = segments.at(seg_index);
          const FT y = seg_data.m_reference_coordinates.y() + qp_result[seg_index];  // differs!
          
          // Check if the angle that seems to be associated to this group of segments is not too close to another value.
          int g_j = -1;
          for (const auto & it_m : values) {
            if (CGAL::abs(it_m.second - y) < m_y_eps) 
              g_j = it_m.first;          
          }
          if (g_j == -1) 
            values[g_i] = y;
          else {                       
            // Merge groups.
            for (const auto gr : groups[g_i]) {
              segments_to_groups_hashmap[gr] = g_j;
              groups[g_j].push_back(gr);
            }
            groups[g_i].clear();
          }
        }
      }

    }

    void assign_segments_to_groups(const std::map <std::size_t, Segment_data> & segments,
                                   std::map<std::size_t, std::vector<std::size_t>> & groups,
                                   std::map<std::size_t, int> & segments_to_groups_hashmap,
                                   std::map<int, FT> & values) {

      for (const auto & sm_i : segments_to_groups_hashmap) {
        int g_i = sm_i.second;
        if (g_i == -1) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = segments.at(seg_index);
          const FT y = seg_data.m_reference_coordinates.y();  // differs!
          int g_j = -1;
          for (const auto & it_m : values) {
            const FT y_j = it_m.second;
            if (CGAL::abs(y_j - y) < m_y_eps) {  // differs!
              g_j = it_m.first;
            }
            if (g_j != -1) 
              break;
          }
          if (g_j == -1) {   
            if (values.size() > 0)             
              g_i = values.rbegin()->first + 1;
            else
              g_i = 0;
            values[g_i] = y;
          } 
          else 
            g_i = g_j;

          segments_to_groups_hashmap[seg_index] = g_i; 
          groups[g_i].push_back(seg_index);
        }
      }

    }

    void build_groups_by_value(const std::map <std::size_t, int> & segments_to_groups_hashmap,
                               const std::map <int, FT> & values, 
                               std::map <FT, std::vector<std::size_t>> & groups_by_value) {

      for (const auto & it_m : values) {
        const FT val = it_m.second;
        if (groups_by_value.find(val) == groups_by_value.end()) 
          groups_by_value[val] = std::vector<std::size_t>();
      }

      for (const auto & sm_i : segments_to_groups_hashmap) {
        const FT val = values.at(sm_i.second);
        if (groups_by_value.find(val) != groups_by_value.end()) 
          groups_by_value[val].push_back(sm_i.first);
      }

    }

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_ORDINATES_2