#ifndef CGAL_VDOL_MISC_H
#define CGAL_VDOL_MISC_H

// this is basically a bunsh of useless functions 
// that should be removed/replaced again in the long term. 

namespace CGAL{
namespace VDOL_3{

template <class Point_3>
typename Point_3::R::FT slength(const Point_3& p){
  return p.x()*p.x()+p.y()*p.y()+p.z()*p.z();
}

template <class Point_3_A, class Point_3_B>
struct Convert_point{
  typedef Point_3_B result_type; 
  Point_3_B operator()(const Point_3_A& a) const {
    return Point_3_B(CGAL::to_double(a.x()),CGAL::to_double(a.y()),CGAL::to_double(a.z()));
  }
};

} // namespace VDOL_3
} // namespace CGAL 

#endif // CGAL_VDOL_MISC_H
