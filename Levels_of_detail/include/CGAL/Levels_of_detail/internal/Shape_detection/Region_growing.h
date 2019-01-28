#ifndef CGAL_LEVELS_OF_DETAIL_REGION_GROWING_H
#define CGAL_LEVELS_OF_DETAIL_REGION_GROWING_H

// STL includes.
#include <queue>
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

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

    using Indices = std::vector<std::size_t>;
    using Regions = std::vector<Indices>;
    
    Region_growing(
      const Indices& indices, 
      Connectivity& connectivity, 
      Conditions& conditions) :
    m_indices(indices),
    m_connectivity(connectivity),
    m_conditions(conditions) { 

      CGAL_precondition(m_indices.size() > 0);
    }

    void detect(Regions &regions) {

      clear();
      regions.clear();
      
      Indices region;
      for (std::size_t i = 0; i < m_indices.size(); ++i) {
        const std::size_t seed_index = m_indices[i];

        // Try to grow a new region from the index of the seed item.
        if (!m_visited[seed_index]) {
          propagate(seed_index, region);

          // Check global conditions.
          if (!m_conditions.is_valid_region(region)) 
            revert(region);
          else 
            regions.push_back(region);
        }
      }
    }

    void clear() {
                
      m_visited.clear();
      m_visited.resize(m_indices.size(), false);
    }

  private:

    void propagate(const std::size_t seed_index, Indices& region) {
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
        Indices neighbors;
        m_connectivity.get_neighbors(item_index, neighbors);

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

    void revert(const Indices& region) {
      for (std::size_t i = 0; i < region.size(); ++i)
        m_visited[region[i]] = false;
    }

    // Fields.
    const Indices& m_indices;
    Input_connectivity& m_connectivity;
    Input_conditions& m_conditions;

    Visited_items m_visited;
    
  }; // Region_growing

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_REGION_GROWING_H
