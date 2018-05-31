#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;

template <typename K, typename Point, typename T>
class Func{
private:
  T* m_tr;
public:
  Func(T* t): m_tr(t){}
  ~Func(){}

  FT operator()(Point query) const{
    return m_tr->compute_func_value_BB(query);
  }
};
