#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <xsimd/xsimd.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Ray_3 Ray_3;
typedef CGAL::Bbox_3  Bbox_3;


struct xRay {
  xsimd::batch<double, 4> originx, originy, originz;
  xsimd::batch<double, 4> directionx, directiony, directionz;
  xsimd::batch<double, 4> inversex;
  xsimd::batch<double, 4> signx;

  xRay(const Ray_3& r)
  {
    originx.broadcast(r.source().x());
    originy.broadcast(r.source().y());
    originz.broadcast(r.source().z());
    directionx.broadcast(r.direction().dx());
    directiony.broadcast(r.direction().dy());
    directionz.broadcast(r.direction().dz());

    inversex.broadcast(1.0);
    inversex  /= directionx;
    // signx =
  }

};


struct xNode {
  xsimd::batch<double, 4> bbxmin, bbxmax, bbymin, bbymax, bbzmin, bbzmax;

  xNode(const Bbox_3 bb0, const Bbox_3 bb1, const Bbox_3 bb2, const Bbox_3 bb3)
  {
    double coord[4] = { bb0.xmin(), bb1.xmin(), bb2.xmin(), bb3.xmin() };
    xsimd::load_aligned(coord, bbxmin);
  }
};


int main(int argc, char* argv[])
{
  Point_3 p(3.1, 3.2, 3.3);
  Vector_3 v(10,11,12);

  Ray_3 r(p,v);

  xRay xr(r);

  Bbox_3 bb = p.bbox();

  xNode xnode(bb,bb,bb,bb);
  return 0;
}
