
// CURRENTLY THIS IS JUST RENAME of the VDOL class. 
// IT SHOULD TAKE A VDOL AND PROVIDE THE CURRENT INTERFACE 
// NEEDED TO GENERATE A MESH 



#ifndef CGAL_VDOL_TO_LABELED_FUNCTION_WRAPPER_3_H
#define CGAL_VDOL_TO_LABELED_FUNCTION_WRAPPER_3_H

#include <CGAL/basic.h>
#include <CGAL/VDOL_3/basic.h> 
#include <CGAL/Voronoi_cell_of_line_3.h>
#include <CGAL/VDOL_3/project_back.h>
#include <boost/variant.hpp>
#include <CGAL/squared_distance_3.h>
#include <CGAL/VDOL_3/Approximation_info.h>


namespace CGAL {


template <class VoronoiDiagramOfLines_3, class LinearKernel>
class VDOL_to_labeled_function_wrapper_3 {
public:
  typedef VoronoiDiagramOfLines_3                    Voronoi_diagram_of_lines_3; 
  typedef LinearKernel                               Linear_kernel; 
  typedef VDOL_to_labeled_function_wrapper_3<Voronoi_diagram_of_lines_3,Linear_kernel> Self; 

public:
  typedef typename Linear_kernel::Point_3   Point_3;
  typedef typename Linear_kernel::Line_3    Line_3;
  typedef typename Linear_kernel::Ray_3     Ray_3;
  typedef typename Linear_kernel::Segment_3 Segment_3; 
  typedef typename Linear_kernel::Sphere_3  Sphere_3;
  typedef typename Linear_kernel::Vector_3  Vector_3;

private:
  typedef typename Voronoi_diagram_of_lines_3::Point_3  VDPoint_3;
  typedef typename Voronoi_diagram_of_lines_3::Line_3   VDLine_3;

private:
  Voronoi_diagram_of_lines_3 m_vdol_3; 
public:
  const Voronoi_diagram_of_lines_3& vdol_3() const { return m_vdol_3;}
  Voronoi_diagram_of_lines_3&       vdol_3()       { return m_vdol_3;}

private:
  VDOL_3::Approximation_info m_approximation_info;
public:
  const VDOL_3::Approximation_info& approximation_info() const {return m_approximation_info;}
  VDOL_3::Approximation_info& approximation_info()             {return m_approximation_info;}
  

public:
//  VDOL_to_labeled_function_wrapper_3(const Voronoi_diagram_of_lines_3& vdol_3) 
//    : m_vdol_3(vdol_3)
//  {}
  
  template <class InputIterator>
  VDOL_to_labeled_function_wrapper_3(InputIterator begin, InputIterator end, int vd_seed = 0)
    :m_vdol_3(begin,end,vd_seed)
  {}
  
  template<class OutputIterator>
  OutputIterator locate(const Point_3& p, OutputIterator oit, int hint = 0) const {
    return this->locate(VDPoint_3(p.x(),p.y(),p.z()), oit, hint);
  }

  double slength(const Point_3& p) const {return p.x()*p.x()+ p.y()*p.y()+p.z()*p.z(); }
  Sphere_3 bounding_sphere() const {return Sphere_3(Point_3(0,0,0),approximation_info().sradius_bounding_sphere);}
  Sphere_3 clipping_sphere() const {return Sphere_3(Point_3(0,0,0),approximation_info().sradius_clipping_sphere);}

  
  int label(const Point_3& p) const {

    if (slength(p)>approximation_info().sradius_clipping_sphere) return 0;

    VDPoint_3 bp(p.x(),p.y(),p.z());
    typedef typename Voronoi_diagram_of_lines_3::FT FT; 

#if 0
    std::vector<int> indices;
    locate(p,std::back_inserter(indices));
    int i = indices[0];
    
    return i+1 + 200 * this->vdol_3().cells()[i].compute_segment_index(bp) ;
#else
    FT min =  CGAL::squared_distance(this->vdol_3().lines()[0],bp);
    int result = 0;
    for(int i = 1;i<this->vdol_3().lines().size();i++){
      FT d = CGAL::squared_distance(this->vdol_3().lines()[i],bp);
      if ( d <= min){
        min = d;
        result = i;
      }
    }
    return result +1 + 200 * this->vdol_3().cells()[result].compute_segment_index(bp) ;
#endif
  }
  
   // Concept Function label 
  typedef int return_type;
  return_type operator()(typename Linear_kernel::Point_3 p, bool b = true ) const {
    return label(p); 
  }    

  template <class OutputIterator> 
  OutputIterator generate_points(OutputIterator oit){

    std::vector<VDPoint_3> vdpoints; 
    vdol_3().generate_points(std::back_inserter(vdpoints),this->approximation_info());

    std::cout << "GEN POINTS SIZE: " << vdpoints.size() << std::endl; 

    std::vector<Point_3> points;     
    std::copy(
        boost::make_transform_iterator(vdpoints.begin(),VDOL_3::Convert_point<VDPoint_3,Point_3>()),
        boost::make_transform_iterator(vdpoints.end()  ,VDOL_3::Convert_point<VDPoint_3,Point_3>()),
        std::back_inserter(points));
    
   
    std::vector<Point_3> filtered_points;   
    for(int i = 0; i < points.size(); i++){
      if( slength(points[i]) < approximation_info().sradius_clipping_sphere + 0.01){
        bool insert = true; 
        for(int j = 0; j < filtered_points.size() && insert ; j++){
          if(CGAL::to_double(CGAL::squared_distance(points[i],filtered_points[j])) < approximation_info().sdistance_final) {
            insert = false ;
          }          
        }
        if(insert) filtered_points.push_back(points[i]);
      }
    }

    std::cout << "Filtered POINTS SIZE: " << filtered_points.size() << std::endl; 
    
    //std::sort(points.begin(),points.end());
    // TODO: remove redundant points 
    // assert(false);
    
    return std::copy(filtered_points.begin(),filtered_points.end(),oit);     
  }
};













} // namespace CGAL 
#endif // CGAL_VDOL_TO_LABELED_FUNCTION_WRAPPER_3_H 
