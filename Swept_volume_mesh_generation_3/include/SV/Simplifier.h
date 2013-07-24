#ifndef SV_SIMPLIFIER_H
#define SV_SIMPLIFIER_H

#include <CGAL/basic.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h> 
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>

namespace SV{

template <class Volume_, class ECM_>
class Simplifier{
public:
  typedef Volume_ Volume;
  typedef ECM_ ECM; 
  
  typedef CGAL::Surface_mesh_simplification::Edge_profile<ECM> Edge_profile; 
  typedef CGAL::Surface_mesh_simplification::Edge_length_cost<ECM> Edge_length_cost;
  typedef CGAL::Surface_mesh_simplification::Midpoint_placement<ECM> Midpoint_placement;
  typedef CGAL::Surface_mesh_simplification::LindstromTurk_placement<ECM> LindstromTurk_placement;
  typedef CGAL::Surface_mesh_simplification::LindstromTurk_cost<ECM> LindstromTurk_cost;
  
  typedef typename Edge_profile::vertex_descriptor         vertex_descriptor;
  typedef typename Edge_profile::edge_descriptor           edge_descriptor;
  typedef typename CGAL::halfedge_graph_traits<ECM>::Point Point; 
  typedef typename Volume::Kernel::FT FT; 
  typedef unsigned int size_type; 
  
  typedef Simplifier<Volume,ECM> Self; 
  
#define SV_SIMPLIFICATION_COST_BOUND 10000.0
private:
  const Volume* m_volume;
public:
  const Volume* volume() const {return m_volume;}
  
public:
  Simplifier(const Volume* volume):m_volume(volume){}
  
public:
  
  class Stop_predicate{
    const Self* m_simplifier; 
  public:
    typedef typename Self::ECM ECM; 
    typedef typename Self::FT  FT; 
    typedef typename Self::size_type size_type; 
    typedef typename Self::Edge_profile Profile; 
    
  public:
    Stop_predicate(const Self* simplifier):m_simplifier(simplifier){}
    const Self* simplifier() const {return m_simplifier;} 
    Self* simplifier() {return m_simplifier;}   
    
    bool operator() (
        const  FT& current_cost,
        const Profile&  profile,
        size_type initial_count,
        size_type current_count) const {
      return current_cost >= SV_SIMPLIFICATION_COST_BOUND;
    }
    
  }; 
  
  
  
  // This is the most important guy
  // If it is not allowed to collapse the edge we have to return
  // an uninitialized boost::optional<FT>
  // otherwise we have to think of a cost function edges are processed 
  // by increasing order 
  class Get_cost{
    const Self* m_simplifier; 
  public:
    typedef typename Self::ECM ECM; 
    typedef typename Self::FT  FT; 
    typedef typename Self::Point Point; 
    typedef typename Self::Edge_profile Profile; 
    typedef boost::optional<FT>  result_type;

  public:
    Get_cost(const Self* simplifier):m_simplifier(simplifier){}
    const Self* simplifier() const {return m_simplifier;} 
    Self* simplifier() {return m_simplifier;}   

    result_type operator()( 
        Profile const& ep, 
        boost::optional<Point> const& placement) const{
      if(CGAL::squared_distance(ep.p0(),ep.p1())>0.1)
        return SV_SIMPLIFICATION_COST_BOUND;

      // we check that the start of the edge is still valid 
      // for the placed vertex. 
      // if so, we return the lengh of the edge
      // otherwise, we return an empty boost::optional 
      
      std::vector<vertex_descriptor> link(ep.link());
      for(int i =0; i<link.size();i++){
        int j=(i+1)%link.size();
        Point p_0=ep.p0();
        Point p_i=link[i]->point();
        Point p_j=link[j]->point();
        
        // compute_facet_badness returning the default constructed 
        // boost optional means that the facet is ok.
        if(simplifier()->volume()->compute_facet_badness(p_0,p_i,p_j))
          return SV_SIMPLIFICATION_COST_BOUND; // result_type();
      }
            
      //return Edge_length_cost()(ep,placement);
      return LindstromTurk_cost()(ep,placement);

      //return FT(1)/CGAL::squared_distance(ep.p0(),ep.p1());
    }
  }; 
  
  // This guy is less important, since it places a new point 
  // we just keep on of the vertices 
  // this way we can use the points to start a new mesh later 
  // we do so by just returning an uninitialized boost::optional 
  class Get_placement{
    const Self* m_simplifier; 
  public:
    typedef typename Self::Point Point; 
    typedef typename Self::Edge_profile Profile; 
    typedef boost::optional<Point>  result_type;
    
    Get_placement(const Self* simplifier):m_simplifier(simplifier){}
    const Self* simplifier() const {return m_simplifier;} 
    Self* simplifier() {return m_simplifier;}   
    
    result_type operator() ( Profile const& edge_profile) const{
      return LindstromTurk_placement()(edge_profile);
      
#if 0 // pull the edge into the vertex with the longest adjacent edge 
      edge_descriptor h0 = edge_profile.v0_v1();
      edge_descriptor h1 = edge_profile.v1_v0();
      
      FT m0(0);  
      edge_descriptor h = h0->next();
      while(h != h1){ 
        FT m0 = m0 + CGAL::squared_distance(h->vertex()->point(),h->opposite()->vertex()->point());
        h = h->opposite()->next();
      }
      
      FT m1(0);  
      h = h1->next();
      while(h != h0){ 
        FT m1 = m1 + CGAL::squared_distance(h->vertex()->point(),h->opposite()->vertex()->point());
        h = h->opposite()->next();
      }
      if(m0<m1){
        return h0->vertex()->point();
      }else{
        return h1->vertex()->point();
      }
#endif

      ///////////////
      return boost::optional<Point>(edge_profile.p0()); // we keep the v_0
      //return Midpoint_placement()(edge_profile); 
      
    }
  };
  
  Get_placement get_placement_object() const {return Get_placement(this);}
  Get_cost get_cost_object() const {return Get_cost(this);}
  Stop_predicate stop_predicate_object() const {return Stop_predicate(this);}
  
};

} // namespace SV 
#endif // SV_SIMPLIFIER_H
