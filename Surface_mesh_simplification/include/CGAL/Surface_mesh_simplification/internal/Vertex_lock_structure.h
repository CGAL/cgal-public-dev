#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/boost/graph/helpers.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tbb/atomic.h>
#include <tbb/concurrent_unordered_map.h>

template<class TM>
class Vertex_lock_structure {
  public:
    typedef boost::graph_traits<TM>                                         GraphTraits;
    typedef boost::graph_traits<const TM>                                   ConstGraphTraits;

    typedef typename GraphTraits::vertex_descriptor                         vertex_descriptor;
    typedef typename GraphTraits::vertex_iterator                           vertex_iterator;
    typedef typename GraphTraits::halfedge_descriptor                       halfedge_descriptor;
    typedef typename GraphTraits::halfedge_iterator                         halfedge_iterator;
    typedef CGAL::Halfedge_around_source_iterator<TM>                       out_edge_iterator;
    typedef CGAL::Halfedge_around_target_iterator<TM>                       in_edge_iterator;
    typedef typename GraphTraits::traversal_category                        traversal_category;
    typedef typename GraphTraits::edges_size_type                           size_type;

    Vertex_lock_structure(const TM& tm) {
        std::pair<vertex_iterator, vertex_iterator> its = vertices(tm);
        while(its.first != its.second) {
            m_locks[*its.first] = false;
            its.first++;
        }
    }

    std::unordered_set<vertex_descriptor>& get_local_locked_set() {
        return m_locked_cells.local();
    }

    bool is_locked_by_this_thread(const vertex_descriptor& o) {
        std::unordered_set<vertex_descriptor>& loc = get_local_locked_set();
        return loc.find(o) != loc.end();
    }

    void set_locked_by_this_thread(const vertex_descriptor& o) {
        std::unordered_set<vertex_descriptor>& loc = get_local_locked_set();
        loc.insert(o);
    }

    bool is_locked(const vertex_descriptor& o) {
        return (m_locks.at(o) == true);
    }

    bool try_lock(const vertex_descriptor& o) {
        return is_locked_by_this_thread(o) 
               || try_lock_impl(o);
    }

    bool try_lock_impl(const vertex_descriptor& o) {
        bool old_value = m_locks.at(o).compare_and_swap(true, false);
        if(old_value == false) {
            set_locked_by_this_thread(o);
            return true;
        }
        return false;
    }

    void unlock_everything_locked_by_this_thread() {
        std::unordered_set<vertex_descriptor>& loc = get_local_locked_set();
        for(const vertex_descriptor& o: loc) {
            m_locks.at(o) = false;
        }
        loc.clear();
    }

    void unlock_all_points_locked_by_this_thread() {
        unlock_everything_locked_by_this_thread();
    }

  protected:
    std::unordered_map<vertex_descriptor, tbb::atomic<bool> > m_locks;  
    tbb::enumerable_thread_specific<std::unordered_set<vertex_descriptor> > m_locked_cells;

};