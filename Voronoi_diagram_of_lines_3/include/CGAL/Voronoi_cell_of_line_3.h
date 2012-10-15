#ifndef CGAL_VORONOI_CELL_OF_LINE_3_H
#define CGAL_VORONOI_CELL_OF_LINE_3_H 

#include <CGAL/VDOL_3/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/envelope_3.h>
#include <CGAL/print_envelop_diagram.h>
#include <CGAL/construct_point_2.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Get_arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>

#include <CGAL/VDOL_3/Single_voronoi_cell_envelope_traits_3.h>
#include <CGAL/Voronoi_cell_of_line_3.h>


#include <CGAL/VDOL_3/vol_tools.h>
#include <CGAL/VDOL_3/svcet_3_state_dependent_functions.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h> 
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <CGAL/VDOL_3/project_back.h> 

namespace CGAL {


template <class LinearKernel, class ArithmeticKernel = CGAL::Arithmetic_kernel>
class Voronoi_cell_of_line_3{

public:
  typedef LinearKernel                                                Linear_kernel; 
  typedef ArithmeticKernel                                            Arithmetic_kernel; 
  typedef Voronoi_cell_of_line_3<Linear_kernel,Arithmetic_kernel>     Self; 
  
private:
  typedef typename Linear_kernel::Line_3    Line_3;
  typedef typename Linear_kernel::Point_3   Point_3; 
  typedef typename Linear_kernel::Vector_3  Vector_3; 
  
  typedef typename Arithmetic_kernel::Integer           Integer; 
  typedef typename Arithmetic_kernel::Rational          Rational; 
  typedef typename Arithmetic_kernel::Bigfloat_interval BFI; 
  
  typedef CGAL::Algebraic_kernel_d_1<Integer>                  AK_1;
  typedef CGAL::Algebraic_curve_kernel_2<AK_1>                 AK_2;
  typedef CGAL::Curved_kernel_via_analysis_2< AK_2 >           CK_2;
  
public: 
  typedef CGAL::VDOL_3::Single_voronoi_cell_envelope_traits_3<Linear_kernel, CK_2>   SVCET_3;

  typedef typename SVCET_3::Poly_int_3 Poly_int_3; 
  typedef typename SVCET_3::Poly_int_2 Poly_int_2; 
  typedef typename SVCET_3::Poly_int_1 Poly_int_1; 

  typedef typename SVCET_3::Poly_rat_3 Poly_rat_3; 
  typedef typename SVCET_3::Poly_rat_2 Poly_rat_2; 
  typedef typename SVCET_3::Poly_rat_1 Poly_rat_1; 
  
  typedef typename SVCET_3::Xy_monotone_surface_3 Xy_monotone_surface_3;
  typedef typename AK_2::Coordinate_1 Coordinate_1 ;

public:
  typedef typename CGAL::Envelope_diagram_2<SVCET_3> Envelope_diagram_2;
private:
  typedef typename Envelope_diagram_2::Topology_traits Topology_traits;
  typedef typename Topology_traits::Vertex Vertex; 
  typedef typename Topology_traits::Halfedge Halfedge;
  typedef typename Topology_traits::Face Face;



  typedef  Arr_naive_point_location<Envelope_diagram_2> Point_location;
  // typedef  Arr_trapezoid_ric_point_location<Envelope_diagram_2> Point_location;
  // typedef  Arr_walk_along_line_point_location<Envelope_diagram_2> Point_location;
  // typedef  Arr_landmarks_point_location<Envelope_diagram_2> Point_location;

class Half_voronoi_cell_of_line_3{
  
private:// member 
  CK_2 m_ckernel; 
  
  std::vector<Coordinate_1>    m_intersection_values; 
  std::vector<Coordinate_1>    m_top_u_values;
  std::vector<Face*>           m_top_faces;
  std::vector<Halfedge*>       m_top_edges;
  
  SVCET_3 m_svcet_3; 
  mutable boost::shared_ptr<Envelope_diagram_2> m_diag; 
  mutable boost::shared_ptr<Point_location> m_pl_diag;


public: // meber access 
  const CK_2& ckernel() const {return m_ckernel;}
  CK_2& ckernel() {return m_ckernel;}

  Linear_kernel ekernel() const {return Linear_kernel();} // TODO

  const SVCET_3& svcet_3() const { return m_svcet_3;}
  SVCET_3& svcet_3(){return m_svcet_3;}
  
  const Envelope_diagram_2& diag() const { return *m_diag;}
  Envelope_diagram_2& diag(){return *m_diag;}

  long number_of_vertical_asymptotes() const {return m_top_edges.size();}
public:
  
  template <class ForwardIterator>
  Half_voronoi_cell_of_line_3(const SVCET_3& svcet_3_, ForwardIterator begin, ForwardIterator end)
    : m_svcet_3(svcet_3_), m_diag(new Envelope_diagram_2(&svcet_3()))
  {
    init_diagram(begin,end);
    init_point_location();
#ifdef CGAL_VDOL_CELL_SAMPLE_CURVES
    init_points();
#endif
  } 
  
  // this function return true if the point is in a cell associated to this line
  // more over it outputs all lines that are in the less equal distance via oit.
  template<class OutputIterator>
  bool contains(const Point_3& p, OutputIterator oit) const {  

    Point_3 lp = in_local_coordinates(p);
    CGAL_precondition(lp.z()>=0);
   
    std::vector<Xy_monotone_surface_3> surfaces;
    
    // get all candidates: 
    if(lp.z() != 0){
      CGAL_precondition(lp.z()>0);
      typename SVCET_3::Point_2 pp = project_point(p);
      compute_surfaces(*m_pl_diag,pp,std::back_inserter(surfaces));
    }else if(lp.y() != 0){
      CGAL_precondition(lp.z()==0);
      CGAL_precondition(lp.y() >0);
      int i = distance( m_top_u_values.begin(), lower_bound(m_top_u_values.begin(),m_top_u_values.end(),lp.x()));
      if(i<m_top_u_values.size() && m_top_u_values[i] == lp.x()){ // take the edge at v[i]
        std::copy(m_top_edges[i]->surfaces_begin(),m_top_edges[i]->surfaces_end(),std::back_inserter(surfaces));
      }else{ // take the face to the left of v[i] 
        std::copy(m_top_faces[i]->surfaces_begin(),m_top_faces[i]->surfaces_end(),std::back_inserter(surfaces));
      }
    }else{
      CGAL_precondition(lp.y()==0);
      CGAL_precondition(lp.z()==0);
      // intersections are vertical lines and can be seen at top too. 
      int i = distance( m_top_u_values.begin(),lower_bound(m_top_u_values.begin(),m_top_u_values.end(),lp.x()));
      // if the point is on such a line report all surfaces on the edge
      if(i<m_top_u_values.size() && m_top_u_values[i] == lp.x()){
        std::copy(m_top_edges[i]->surfaces_begin(),m_top_edges[i]->surfaces_end(),std::back_inserter(surfaces));
      }
      // otherwise the point is just on the line and we do nothing 
    }
    
    // compute the minimal distance of p to base line and all reported surfaces. 
    Rational min_sdist = CGAL::squared_distance(svcet_3().base_line_3(),p);
    for(int i=0; i<surfaces.size(); i++){
      min_sdist = (std::min)(min_sdist,CGAL::squared_distance(surfaces[i].line(),p));
    }
    
    // report all surfaces that are have the minimal distance 
    for(int i=0; i<surfaces.size(); i++){
      if(CGAL::squared_distance(surfaces[i].line(),p) == min_sdist)
        *oit++ = surfaces[i].line();
    }
    
    // also report the base line if it has the minimal distance 
    if (CGAL::squared_distance(svcet_3().base_line_3(),p) == min_sdist){
      *oit++ = svcet_3().base_line_3();
      return true;
    }
    
    return false; 
  }
 

  // compute the index of the segment the point is associated to 
  // relavant in case the line is intersected by other lines. 
  int compute_segment_index(const Point_3& p) const {
    if(m_intersection_values.size()==0) return 0; 
    Rational  t =  in_local_coordinates(p).x() ;
    return distance(
        m_intersection_values.begin(),
        lower_bound(m_intersection_values.begin(),m_intersection_values.end(),t));
  }
  
  void init_point_location(){
    // point location for normal points
    m_pl_diag = boost::shared_ptr<Point_location>(new Point_location(*m_diag));

    // point location for points in H^*
    Vertex   *tv = m_diag->topology_traits()->top_left_vertex();
    Halfedge *ht = tv->halfedge();
    if(VDOL_3::is_vertical(ht)) ht = ht->next()->opposite();
    CGAL_assertion(!VDOL_3::is_vertical(ht));
    CGAL_assertion(ht->vertex() == tv);
    
    ht = ht->opposite()->next()->opposite();
    while(ht->vertex() != m_diag->topology_traits()->top_right_vertex()){
      CGAL_assertion(ht->vertex() != m_diag->topology_traits()->top_right_vertex());
      CGAL_assertion(ht->vertex() != m_diag->topology_traits()->top_left_vertex());
      CGAL_assertion(ht->vertex() != m_diag->topology_traits()->bottom_right_vertex());
      CGAL_assertion(ht->vertex() != m_diag->topology_traits()->bottom_left_vertex());
      CGAL_assertion(VDOL_3::degree(ht->vertex())==3);

      m_top_edges.push_back(ht->next());
      m_top_faces.push_back(ht->next()->opposite()->outer_ccb()->face());

      // m_top_u_values.push_back(ht->vertex()->point().x()); 
      m_top_u_values.push_back(
          ht->next()->opposite()->direction() == ARR_LEFT_TO_RIGHT? 
          ht->next()->opposite()->curve().curve_end_x(ARR_MAX_END):
          ht->next()->opposite()->curve().curve_end_x(ARR_MIN_END));
      
      // vertical curves are only due to intersections with another line 
      if(ht->next()->curve().is_vertical())
        m_intersection_values.push_back(ht->next()->opposite()->curve().x());
 
      ht = ht->opposite()->next()->opposite();
    } 
    m_top_faces.push_back(ht->outer_ccb()->face());
    std::cerr <<" top_vertices.size() " << m_top_u_values.size() <<std::endl;
    std::cerr <<" intersection_values.size() " << m_intersection_values.size() <<std::endl;
  }
  
  template <class ForwardIterator>
  void init_diagram(ForwardIterator begin, ForwardIterator end){
    std::cerr<< " Start Envelope_3 " << std::endl;
    CGAL::Timer timer_1;
    timer_1.start();
    CGAL::lower_envelope_3(begin, end, *m_diag);
    timer_1.stop();
    std::cerr<< " Time: " << timer_1.time() << std::endl;
    print_envelop_diagram(*m_diag);
  }
  
  Point_3 in_local_coordinates(const Point_3& ep) const {
    Vector_3 ev_diff = ekernel().construct_vector_3_object()(svcet_3().origin(),ep);
    Rational  lv1 = ekernel().compute_squared_length_3_object()(svcet_3().v1());
    Rational  lv2 = ekernel().compute_squared_length_3_object()(svcet_3().v2());
    Rational  lv3 = ekernel().compute_squared_length_3_object()(svcet_3().v3());
    
    Rational  c1 = ekernel().compute_scalar_product_3_object()(svcet_3().v1(),ev_diff) / lv1; 
    Rational  c2 = ekernel().compute_scalar_product_3_object()(svcet_3().v2(),ev_diff) / lv2; 
    Rational  c3 = ekernel().compute_scalar_product_3_object()(svcet_3().v3(),ev_diff) / lv3; 
    return Point_3(c1,c2,c3);
  }

  typename SVCET_3::Point_2 project_point(const Point_3& ep) const {
    Point_3 lep = in_local_coordinates(ep);
    CGAL_precondition(lep.z() != 0);
    Rational u = lep.x(); 
    Rational v = lep.y()/lep.z();
    return CGAL::construct_point_2(&svcet_3() ,Rational(u),Rational(v));
  }


  // APROXIMATION OF VERTICES AND EDGES 
  template <class OutputIterator>
  OutputIterator generate_vertex_points ( OutputIterator oit, const VDOL_3::Approximation_info& info) const {
    typedef typename Envelope_diagram_2::Vertex_const_iterator  VIT;
    for(VIT vit = this->diag().vertices_begin(); vit != this->diag().vertices_end(); vit++){
      oit = project_back(svcet_3(),*vit,info,oit);  
    }
    return oit; 
  }

  template <class OutputIterator>
  OutputIterator generate_edge_points ( OutputIterator oit, const VDOL_3::Approximation_info& info) const {
    typedef typename Envelope_diagram_2::Edge_const_iterator EIT;
    for(EIT eit = this->diag().edges_begin(); eit != this->diag().edges_end();eit++){
      oit = project_back(svcet_3(),*eit,info,oit);  
    }
    return oit; 
  } 

  template <class OutputIterator>
  OutputIterator generate_face_points ( OutputIterator oit, const VDOL_3::Approximation_info& info) const {
    CGAL_SNAP_SVCET_3_TYPEDEFS;
    
    // construct polynomial representing clippling sphere 
    Poly_int_3 x = shift(Poly_int_3(1),1,0);
    Poly_int_3 y = shift(Poly_int_3(1),1,1);
    Poly_int_3 z = shift(Poly_int_3(1),1,2);
    Poly_int_3 csphere = x*x+y*y+z*z - info.sradius_clipping_sphere; 
    Poly_int_3 tcsphere = VDOL_3::construct_transformed_bisector_3(&m_svcet_3,csphere);   

    
    // iterator over all faces and collect all bisector surfaces 
    typedef typename Envelope_diagram_2::Face_iterator FIT;
    //std::cerr << " COLLECT SURFACES "  << std::endl;
    std::set<Xy_monotone_surface_3> surfaces;
    for(FIT fit = m_diag->faces_begin(); fit != m_diag->faces_end();fit++){
      if(fit->number_of_surfaces() > 0 )
        surfaces.insert(fit->surface());
    }
    //std::cerr << " NUMBER OF BISECTORS "  << std::endl;
    
    // intersect bisector surfaces with clipping sphere 
    std::vector<typename SVCET_3::Curve_2> curves;
    for(typename std::set<Xy_monotone_surface_3>::iterator sit = surfaces.begin(); sit != surfaces.end();sit++){
      Poly_int_2 projection = canonicalize(make_square_free(resultant(sit->transformed_bisector(),tcsphere)));
      curves.push_back(svcet_3().kernel().construct_curve_2_object()(projection));
    }
    // split prjected curves into x_monotone non intersecting subcurves 
    std::vector<typename SVCET_3::X_monotone_curve_2> subcurves;
    // Though, SVCET_3 is a valid traits vor Arr_2, it is important to use the curved_kernel 
    // since Intersect_2 of SVCET_3 may throw which we don't want here. 
    compute_subcurves (curves.begin(),curves.end(),std::back_inserter(subcurves),false,svcet_3().curved_kernel_2());
    
    typedef typename std::vector<typename SVCET_3::X_monotone_curve_2>::iterator SCIT;  
    typename Envelope_diagram_2::Face_const_handle   fh;
    for(SCIT scit= subcurves.begin();scit != subcurves.end(); scit++){   
      std::cerr << "."  << std::flush; 
      Point_2 p2 = svcet_3().construct_interior_vertex_2_object()(*scit);
      if(CGAL::assign (fh,m_pl_diag->locate(p2))){
        if(fh->number_of_surfaces() > 0){
          // do not use tcsphere to project back since this is unstable 
          typename Coercion_traits<Poly_int_3,Poly_rat_3>::Cast cast;
          Point_3 B = project_back(svcet_3(),p2,fh->surface().transformed_bisector());
          if(CGAL::abs(CGAL::to_double(canonicalize(cast(csphere)).evaluate(B.z()).evaluate(B.y()).evaluate(B.x()))) < 0.0001){
            project_back(svcet_3(),*scit,fh->surface().transformed_bisector(),info,oit);
          }
        }
      }
    }
    std::cerr << "."  << std::endl;
    return oit; 
  }
  
  // reports total number of vertices (first) and (second) number of vertices
  // that are really relevant for the voronoi diagram. 
  std::pair<int,int> number_of_vertices() const {
    int number_of_needed_vertices(0);
    typedef typename Envelope_diagram_2::Vertex_iterator VIT; 
    for(VIT vit = m_diag->vertices_begin(); vit != m_diag->vertices_end(); vit++){
      if(vit->degree()>=3) number_of_needed_vertices++;
    }
    return std::make_pair(
        m_diag->number_of_vertices(),
        number_of_needed_vertices);
  }
  
  template< class OutputIterator> 
  OutputIterator report_bisectors(OutputIterator oit) const {
    typedef typename Envelope_diagram_2::Face_iterator FIT; 
    typedef typename Envelope_diagram_2::Surface_const_iterator SIT;
    for(FIT fit = m_diag->faces_begin(); fit != m_diag->faces_end();fit++){
      for(SIT sit = fit->surfaces_begin(); sit != fit->surfaces_end();sit++){
        *oit++ = sit->bisector();
      }
    }
    return oit; 
  }

    
  template< class OutputIterator> 
  OutputIterator report_transformed_bisectors(OutputIterator oit) const {
    typedef typename Envelope_diagram_2::Face_iterator FIT; 
    typedef typename Envelope_diagram_2::Surface_const_iterator SIT;
    for(FIT fit = m_diag->faces_begin(); fit != m_diag->faces_end();fit++){
      for(SIT sit = fit->surfaces_begin(); sit != fit->surfaces_end();sit++){
        *oit++ = sit->transformed_bisector();
      }
    }
    return oit; 
  }
  
  template< class OutputIterator> 
  OutputIterator report_resultants(OutputIterator oit) const {
    typedef typename Envelope_diagram_2::Vertex_iterator VIT; 
    typedef typename Envelope_diagram_2::Surface_const_iterator SIT;
    for(VIT vit = m_diag->vertices_begin(); vit != m_diag->vertices_end();vit++){
      Poly_int_1 poly = vit->point().x().polynomial();
      *oit++ = poly;
    }
    return oit; 
  }
  

}; 
  
//  MEMBERS 
  Half_voronoi_cell_of_line_3 m_hvcl_3_pos;
  Half_voronoi_cell_of_line_3 m_hvcl_3_neg;
public: 
 

  template <class ForwardIterator>
  Voronoi_cell_of_line_3(int seed, ForwardIterator begin, ForwardIterator end)
    :m_hvcl_3_pos(SVCET_3(*begin++,seed),begin,end),
     m_hvcl_3_neg(m_hvcl_3_pos.svcet_3().mirror(),begin,end)
  {}

  // use positive side 
  bool compute_side(const Point_3& ep) const {
    Linear_kernel ek =  m_hvcl_3_pos.ekernel();
    SVCET_3 svcet_3 = m_hvcl_3_pos.svcet_3();
    Vector_3 ev_diff =ek.construct_vector_3_object()(svcet_3.origin(),ep);
    Sign      s1  = CGAL::sign(ek.compute_scalar_product_3_object()(svcet_3.v3(),ev_diff));
    Sign      s2  = CGAL::sign(ek.compute_scalar_product_3_object()(svcet_3.v2(),ev_diff));
    if(s1 == CGAL::POSITIVE) return true;
    if(s1 == CGAL::NEGATIVE) return false;
    if(s2 == CGAL::POSITIVE) return true;
    if(s2 == CGAL::NEGATIVE) return false;
    return true;
  }

 
  template<class OutputIterator>
  bool contains(const Point_3& ep, OutputIterator oit) const {
#ifdef CGAL_VDOL_USE_COUNTER
    count_cell_contains++;
#endif
    return compute_side(ep)? 
      m_hvcl_3_pos.contains(ep , oit):
      m_hvcl_3_neg.contains(ep , oit);
  }
  
  int compute_segment_index(const Point_3& p) const {
    return m_hvcl_3_pos.compute_segment_index(p);
  }
  
  const Line_3& line_3() const {return m_hvcl_3_pos.svcet_3().base_line_3();}  

  int number_of_vertical_asymptotes() const { 
    return m_hvcl_3_pos.number_of_vertical_asymptotes();
  }
  
  const Envelope_diagram_2& envelope_diagram_2(bool b) const {
    return (b) ? m_hvcl_3_pos.diag() : m_hvcl_3_neg.diag();
  }


  
  template <class OutputIterator>
  OutputIterator generate_vertex_points(
      OutputIterator oit, const VDOL_3::Approximation_info& approximation_info ) const {
    oit = m_hvcl_3_pos.generate_vertex_points(oit,approximation_info);
    oit = m_hvcl_3_neg.generate_vertex_points(oit,approximation_info);
    return oit; 
  }
  template <class OutputIterator>
  OutputIterator generate_edge_points(
      OutputIterator oit, const VDOL_3::Approximation_info& approximation_info ) const {
    oit = m_hvcl_3_pos.generate_edge_points(oit,approximation_info);
    oit = m_hvcl_3_neg.generate_edge_points(oit,approximation_info);
    return oit; 
  } 
  template <class OutputIterator>
  OutputIterator generate_face_points(
      OutputIterator oit, const VDOL_3::Approximation_info& approximation_info ) const {
    oit = m_hvcl_3_pos.generate_face_points(oit,approximation_info);
    oit = m_hvcl_3_neg.generate_face_points(oit,approximation_info);
    return oit; 
  }
  
  // reports total number of vertices (first) and (second) number of vertices
  // that are really relevant for the voronoi diagram.
  std::pair<int,int> number_of_vertices() const {
    std::pair<int,int> pos = m_hvcl_3_pos.number_of_vertices(); 
    std::pair<int,int> neg = m_hvcl_3_neg.number_of_vertices(); 
    return std::make_pair( pos.first + neg.first , pos.second + neg.second ); 
  }

  template<class OutputIterator> 
  OutputIterator report_bisectors(OutputIterator oit) const {
    oit = m_hvcl_3_pos.report_bisectors(oit);
    oit = m_hvcl_3_neg.report_bisectors(oit);
    return oit; 
  }
  template<class OutputIterator> 
  OutputIterator report_transformed_bisectors(OutputIterator oit) const {
    oit = m_hvcl_3_pos.report_transformed_bisectors(oit);
    oit = m_hvcl_3_neg.report_transformed_bisectors(oit);
    return oit; 
  }
  template<class OutputIterator> 
  OutputIterator report_resultants(OutputIterator oit) const {
    oit = m_hvcl_3_pos.report_resultants(oit);
    oit = m_hvcl_3_neg.report_resultants(oit);
    return oit; 
  }
  
};

} // namespace CGAL
 
#endif // CGAL_VORONOI_CELL_OF_LINE_3_H
