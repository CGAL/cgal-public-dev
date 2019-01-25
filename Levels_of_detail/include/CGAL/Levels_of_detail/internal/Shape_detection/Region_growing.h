#ifndef CGAL_LEVELS_OF_DETAIL_REGION_GROWING_H
#define CGAL_LEVELS_OF_DETAIL_REGION_GROWING_H

// STL includes.
#include <queue>
#include <vector>

// CGAL includes.
#include <CGAL/Iterator_range.h>

namespace CGAL {
namespace Levels_of_detail {

  template<
  typename InputRange, 
  typename Connectivity, 
  typename Conditions>
  class Region_growing {

  public:
    using Input_range = InputRange;
    using Input_connectivity = Connectivity;
    using Input_conditions = Conditions;

    using Visited_items = std::vector<bool>;
    using Running_queue = std::queue<std::size_t>;    

    using Items = std::vector<std::size_t>;
    using Regions = std::vector<Items>;
    using Region_range = CGAL::Iterator_range<typename Regions::const_iterator>;
    using Item_range = CGAL::Iterator_range<typename Items::const_iterator>;
    
    Region_growing(
      const InputRange& input_range, 
      Connectivity& connectivity, 
      Conditions& conditions) :
    m_input_range(input_range),
    m_connectivity(connectivity),
    m_conditions(conditions)
    { }

    void detect() {

      clear();
      Items region;

      for (std::size_t seed_index = 0; 
        seed_index < m_input_range.size(); 
        ++seed_index) {
                    
        // Try to grow a new region from the index of the seed item.
        if (!m_visited[seed_index]) {
          propagate(seed_index, region);

          // Check global conditions.
          if (!m_conditions.is_valid_region(region)) 
            revert(region);
          else 
            m_regions.push_back(region);
        }
      }
      m_output_regions = 
      Region_range(m_regions.begin(), m_regions.end());

      // Return indices of all unassigned items.
      for (std::size_t item_index = 0; 
        item_index < m_input_range.size(); 
        ++item_index) {
                    
        if (!m_visited[item_index]) 
          m_unassigned.push_back(item_index);
      }

      m_output_unassigned 
      = Item_range(m_unassigned.begin(), m_unassigned.end());
    }

    const Region_range& regions() const {
      return m_output_regions;
    }

    const Item_range& unassigned_items() const {
      return m_output_unassigned;
    }

    const std::size_t number_of_regions() const {
      return m_regions.size();
    }

    const std::size_t number_of_unassigned_items() const {
      return m_unassigned.size();
    }

    void clear() {
                
      m_visited.clear();
      m_regions.clear();

      m_unassigned.clear();
      m_visited.resize(m_input_range.size(), false);
    }

  private:

    void propagate(const std::size_t seed_index, Items& region) {
      region.clear();

      // Use two queues, while running on this queue, push to the other queue;
      // When the queue is done, update the shape of the current region and swap to the other queue;
      // depth_index is the index of the queue we are using.
      Running_queue running_queue[2];
      bool depth_index = 0;

      // Once the index of an item is pushed to the queue, it is pushed to the region too.
      m_visited[seed_index] = true;
      running_queue[depth_index].push(seed_index);
      region.push_back(seed_index);

      // Update internal properties of the propagating region.
      m_conditions.update(region);

      while (
        !running_queue[depth_index].empty() || 
        !running_queue[!depth_index].empty()) {

        // Call the next item index of the queue and remove it from the queue.
        const std::size_t item_index = running_queue[depth_index].front();
        running_queue[depth_index].pop();

        // Get neighbors of the current item.
        Items neighbors;
        m_connectivity.get_neighbors(item_index, std::back_inserter(neighbors));

        // Visit the neighbors.
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
          const std::size_t neighbor_index = neighbors[i];
                        
          if (
            !m_visited[neighbor_index] && 
            m_conditions.belongs_to_region(neighbor_index, region)) {

            // Add this neighbor to the other queue so that we can visit it later.
            m_visited[neighbor_index] = true;
            running_queue[!depth_index].push(neighbor_index);
            region.push_back(neighbor_index);
          }
        }

        // Update internal properties of the propagating region.
        if (running_queue[depth_index].empty()) {

          m_conditions.update(region);
          depth_index = !depth_index;
        }
      }
    }

    void revert(const Items& region) {
      for (std::size_t i = 0; i < region.size(); ++i)
        m_visited[region[i]] = false;
    }

    // Fields.
    const Input_range& m_input_range;
    Input_connectivity& m_connectivity;
    Input_conditions& m_conditions;

    Visited_items m_visited;
    Regions m_regions;
    Items m_unassigned;

    Region_range m_output_regions = 
    Region_range(m_regions.begin(), m_regions.end());

    Item_range m_output_unassigned = 
    Item_range(m_unassigned.begin(), m_unassigned.end());
  };

} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_REGION_GROWING_H
