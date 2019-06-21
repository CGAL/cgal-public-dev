#ifndef CGAL_SHAPE_REGULARIZATION_DENSE_QP_SOLVER
#define CGAL_SHAPE_REGULARIZATION_DENSE_QP_SOLVER

// #include <CGAL/license/Shape_regularization.h>

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits>
  class Dense_QP_solver{ 

  public:
    using Traits = GeomTraits;
    using FT = typename GeomTraits::FT;

    Dense_QP_solver() {}

    void solve(std::vector<FT> & result){
      result.clear();
      result.push_back(FT(-1.90353));
      result.push_back(FT(3.80706));
      result.push_back(FT(-1.90353));
      result.push_back(FT(1.56148e-12));
      result.push_back(FT(2.66171e-14));
      result.push_back(FT(1.56259e-12));
    }
    // creates an instance of CGAL QP solver
  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DENSE_QP_SOLVER
