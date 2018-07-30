

#include <CGAL/Simple_cartesian.h>
typedef double FT;
typedef CGAL::Simple_cartesian<FT> K;

typedef CGAL::Simple_cartesian<CGAL::Interval_nt_advanced> FK;
typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;
typedef CGAL::Cartesian_converter<K, EK> C2E;
typedef CGAL::Cartesian_converter<K, FK> C2F;

typedef CGAL::Filtered_predicate<EK::Orientation_2,
                                 FK::Orientation_2, C2E, C2F> Filtered_Orientation_2;

typedef K::Point_2 Point;
typedef EK::Point_2 ePoint;

using boost::math::float_advance;

CGAL::Orientation
exact_orientation(const Point& , const Point&, const Point& r)
{
  if(r.x() == r.y()) return CGAL::COLLINEAR;
  if(r.x() > r.y()) return CGAL::LEFT_TURN;
  return CGAL::RIGHT_TURN;
}

int main()
{
  std::cout.precision(17);
  Point p(0.5, 0.5), q(24.,24.);

  FT ry = 0.5;
  for(int i = 0; i < 100; i++){
    FT rx = 0.5;
    for(int j = 0; j < 100; j++){
      Point r(rx,ry);
      CGAL::Orientation eori = Filtered_Orientation_2()(q,p,r);
      CGAL::Orientation eori2 = exact_orientation(q,p,r);
      if(eori != eori2){
        std::cout << eori << " != " << eori2 << std::endl;
      }

#if 0
      CGAL::Orientation ori = K::Orientation_2()(q,p,r);
      if(ori == eori){
        std::cout << " ";
      }else{
        std::cout << "x";
      }
#endif      
      rx = float_advance(rx,1);
    }
    ry = float_advance(ry,1);

  }
  return 0;
}
