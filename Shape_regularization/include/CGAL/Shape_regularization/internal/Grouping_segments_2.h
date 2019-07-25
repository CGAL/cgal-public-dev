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
    typename GeomTraits>
  class Grouping_segments_2 {
  public:
    using Traits = GeomTraits;
    using FT = typename GeomTraits::FT;
    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Targets_map = std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>>;
    using Relations_map = std::map <std::pair<std::size_t, std::size_t>, std::pair<int, std::size_t>>;
    enum Regularization_type {ANGLES, ORDINATES};

    Grouping_segments_2(
      Regularization_type type) :
    m_type(type),
    m_tolerance(FT(1) / FT(1000000)) {

      CGAL_precondition(m_type == ANGLES || m_type == ORDINATES);

      switch (m_type) {
        case ANGLES:
          m_eps = FT(1) / FT(4);
          break;
        case ORDINATES:
          m_eps = FT(1);
          break;
      }

    }

      void make_groups(const std::size_t n, const std::map <std::size_t, Segment_data> & segments,
                       const std::vector<FT> & qp_result,
                       std::map<FT, std::vector<std::size_t>> & groups_by_value,
                       const Targets_map & t_ijs, const Relations_map & r_ijs = Relations_map()) { 

      CGAL_precondition(n > 0);
      CGAL_precondition(qp_result.size() > 0);
      groups_by_value.clear();
      
      m_groups.clear();
      m_segments_to_groups_hashmap.clear();

      for (const auto & it : segments) {
        std::size_t seg_index = it.second.m_index;
        m_segments_to_groups_hashmap[seg_index] = -1;
      }

      build_initial_groups(n, t_ijs, r_ijs, segments, qp_result);

      std::map<int, FT> values;

      build_map_of_values(qp_result, segments, values);

      // Try to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group.
      assign_segments_to_groups(segments, values);

      build_groups_by_value(values, groups_by_value);

    } 

  private:
    FT m_eps;
    const FT m_tolerance;
    Regularization_type m_type;
    std::map<std::size_t, int> m_segments_to_groups_hashmap;
    std::map <std::size_t, std::vector<std::size_t>> m_groups;
    
    void build_initial_groups(const std::size_t n,
                              const std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>> & t_ijs,
                              const std::map <std::pair<std::size_t, std::size_t>, std::pair<int, std::size_t>> & r_ijs,
                              const std::map <std::size_t, Segment_data> & segments,
                              const std::vector<FT> & qp_result) {

      std::size_t g = 0;
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
        CGAL_postcondition(r == 0 || r == 1);

        if (CGAL::abs(qp_result[n + p]) >= m_tolerance) { 
          if(rel_it != r_ijs.end()) ++rel_it;
          continue;
        }

        const int g_i = m_segments_to_groups_hashmap[i];
        const int g_j = m_segments_to_groups_hashmap[j];
        const int groups_status = check_group_status(g_i, g_j);

        switch (groups_status) {
          case -1: break;

          case 1:
            r == 0 ? create_single_group(i, j, g) : create_separate_groups(i, j, g);
            break;

          case 2:
            r == 0 ? assign_segment_to_group(i, j) : create_new_group(i, g);
            break;

          case 3:
            r == 0 ? assign_segment_to_group(j, i) : create_new_group(j, g);
            break;

          case 4:
            if (r == 0) merge_two_groups(g_i, g_j);
            break;
        }
  
        if(rel_it != r_ijs.end()) ++rel_it;
      }

    }

    void build_map_of_values(const std::vector<FT> & qp_result,
                             const std::map <std::size_t, Segment_data> & segments,
                             std::map<int, FT> & values) {

      for (const auto & sm_i : m_segments_to_groups_hashmap) {
        int g_i = sm_i.second;

        if (g_i != -1 && (values.find(g_i) == values.end())) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = segments.at(seg_index);
          const FT val = get_value_1(seg_data, qp_result);
          
          // Check if the angle that seems to be associated to this group of segments is not too close to another value.
          int g_j = -1;
          for (const auto & it_m : values) {
            if (CGAL::abs(it_m.second - val) < m_eps) 
              g_j = it_m.first;
          }

          if (g_j == -1) 
            values[g_i] = val;
          else                       
            merge_two_groups(g_j, g_i);

        }
      }

    }

    void assign_segments_to_groups(const std::map <std::size_t, Segment_data> & segments,
                                   std::map<int, FT> & values) {

      for (const auto & sm_i : m_segments_to_groups_hashmap) {
        int g_i = sm_i.second;
        if (g_i == -1) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = segments.at(seg_index);
          const FT val = get_value_2(seg_data);
          int g_j = -1;

          for (const auto & it_m : values) {
            const FT val_j = it_m.second;
            g_j = get_g_j(val, val_j, it_m);
            if (g_j != -1) 
              break;
          }

          if (g_j == -1) {   
            values.size() > 0 ? g_i = values.rbegin()->first + 1 : g_i = 0;
            values[g_i] = val;
          } 
          else 
            g_i = g_j;

          m_segments_to_groups_hashmap[seg_index] = g_i;
          m_groups[g_i].push_back(seg_index);
        }
      }

    }

    void build_groups_by_value(std::map <int, FT> & values, 
                               std::map <FT, std::vector<std::size_t>> & groups_by_value) {

      for (const auto & it_m : values) {
        const FT val = it_m.second;
        if (groups_by_value.find(val) == groups_by_value.end()) 
          groups_by_value[val] = std::vector<std::size_t>();
      }

      for (const auto & sm_i : m_segments_to_groups_hashmap) {
        // If segment s_i is included in a group of parallel segments,
        // then it should be assigned to a leaf of the regularization tree.
        const FT val = values.at(sm_i.second);
        if (groups_by_value.find(val) != groups_by_value.end()) 
          groups_by_value[val].push_back(sm_i.first);
      }

    } 

    int check_group_status(const int g_i, const int g_j) {

      if (g_i == -1 && g_j == -1) return 1;
      if (g_i == -1 && g_j != -1) return 2;
      if (g_i != -1 && g_j == -1) return 3;
      if (g_i != -1 && g_j != -1 && g_i != g_j) return 4;
      return -1;

    }

    void create_single_group (const std::size_t i, const std::size_t j, std::size_t & g) {

      m_segments_to_groups_hashmap[i] = g;
      m_segments_to_groups_hashmap[j] = g;

      m_groups[g].push_back(i);
      m_groups[g].push_back(j); 

      ++g;
    }

    void create_separate_groups(const std::size_t i, const std::size_t j, std::size_t & g) {

      create_new_group(i, g);
      create_new_group(j, g);

    }

    void assign_segment_to_group(const std::size_t i, const std::size_t j) {

      const int g_j = m_segments_to_groups_hashmap[j];
      m_segments_to_groups_hashmap[i] = g_j;
      m_groups[g_j].push_back(i);

    }

    void create_new_group(const std::size_t i, std::size_t & g) {

      m_segments_to_groups_hashmap[i] = g;
      m_groups[g].push_back(i);
      ++g;

    }

    void merge_two_groups(const int g_i, const int g_j) {

      for (const auto gr : m_groups[g_j]) {
        m_segments_to_groups_hashmap[gr] = g_i;
        m_groups[g_i].push_back(gr);
      }
      m_groups[g_j].clear();
      
    }

    FT get_value_1(const Segment_data & seg_data,
                     const std::vector<FT> & qp_result) {

      const std::size_t seg_index = seg_data.m_index;
      FT val = 0;

      switch (m_type) {
        case ANGLES:
        {
          val = seg_data.m_orientation + qp_result[seg_index]; 

          if (val < FT(0)) val += FT(180); 
          else if (val > FT(180)) val -= FT(180);

          break;
        }

        case ORDINATES:
          val = seg_data.m_reference_coordinates.y() + qp_result[seg_index];
          break;
      }
      return val;

    }

     FT get_value_2(const Segment_data & seg_data) {

       FT val = 0;
       switch (m_type) {
        case ANGLES:
          val = seg_data.m_orientation;
          break;

        case ORDINATES:
          val = seg_data.m_reference_coordinates.y();
          break;
      }

      return val;
     }

      int get_g_j(const FT val, const FT val_j, const auto & it_m) {

       int g_j = -1;
       switch (m_type) {
        case ANGLES:
        {
          for (int k = -1; k <= 1; ++k) {  
            if (CGAL::abs(val_j - val + static_cast<FT>(k) * FT(180)) < m_eps) {  
              g_j = it_m.first;
              break;  
            }
          }
          break;
        }

        case ORDINATES:
          if (CGAL::abs(val_j - val) < m_eps)  
            g_j = it_m.first;
          break;
       }

      return g_j;
     }


  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2